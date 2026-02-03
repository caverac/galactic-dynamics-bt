"""Module for multipole expansion of galaxy density models."""

from abc import ABC, abstractmethod
from dataclasses import dataclass
import logging
from pathlib import Path
from typing import Generic, TypeVar

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d

from galactic_dynamics_bt.utils.assets import register_asset, save_figure_if_changed

T = TypeVar("T")
logger = logging.getLogger(__name__)


class Model(ABC, Generic[T]):
    """Abstract base class for galaxy density and potential models."""

    params: T

    def __init__(self, params: T) -> None:
        """Initialize the model with given parameters."""
        self.params = params

    @abstractmethod
    def density(self, R: float, z: float) -> float:
        """Calculate the density at cylindrical coordinates (R, z)."""

    @abstractmethod
    def potential(self, R: float, z: float) -> float:
        """Calculate the gravitational potential at cylindrical coordinates (R, z)."""


@dataclass
class SatohModelParams:
    """Parameters for the Satoh model."""

    q: float


class SatohModel(Model[SatohModelParams]):
    """Satoh (1980) model for flattened galaxy density profiles."""

    params: SatohModelParams

    def density(self, R: float, z: float) -> float:
        """
        Calculate the density at cylindrical coordinates (R, z) using the Satoh model.

        Parameters
        ----------
        R : float
            The cylindrical radius.
        z : float
            The height above the plane.

        Returns
        -------
        float
            The density at the given coordinates.
        """
        q = self.params.q

        l2 = z**2 + q**2
        S2 = R**2 + z**2 + (1 + 2 * np.sqrt(l2))
        S3 = S2**1.5

        rho = (q**2 / (4 * np.pi * S3 * l2)) * (1 / np.sqrt(l2) + 3 * (1 - (R**2 + z**2) / S2))
        return float(rho)

    def potential(self, R: float, z: float) -> float:
        """
        Calculate the gravitational potential at cylindrical coordinates (R, z) using the Satoh model.

        Parameters
        ----------
        R : float
            The cylindrical radius.
        z : float
            The height above the plane.

        Returns
        -------
        float
            The gravitational potential at the given coordinates.
        """
        q = self.params.q

        l2 = z**2 + q**2
        S = np.sqrt(R**2 + z**2 + (1 + 2 * np.sqrt(l2)))

        phi = -1 / S
        return float(phi)


class RhoL0Interpolator:
    """
    Interpolator class for pre-computing rho_l0 values for efficient integration.

    This class pre-computes the l-th multipole moment of the density at various
    radii and creates an interpolator for fast evaluation during integration.
    """

    def __init__(
        self,
        model: Model,
        multipole_order: int,
        /,
        *,
        Ngrid: int = 500,
        a_min: float = 1e-5,
        a_max: float = 1e7,
    ) -> None:
        """
        Initialize the interpolator for a specific multipole order l.

        Parameters
        ----------
        model : Model
            The galaxy density model.
        multipole_order : int
            The multipole order.
        Ngrid : int, optional
            Number of grid points for pre-computation, by default 500.
        """
        self.model = model
        self.multipole_order = multipole_order
        self.Ngrid = Ngrid

        self.a_grid = np.logspace(np.log10(a_min), np.log10(a_max), Ngrid)
        self.rho_l0_values = np.zeros(Ngrid)
        self.a_min = a_min
        self.a_max = a_max

        # Pre-compute rho_l0 values at grid points
        for i, a in enumerate(self.a_grid):
            self.rho_l0_values[i] = self._compute_rho_l0(a)

        self.interpolator = interp1d(
            self.a_grid,
            self.rho_l0_values,
            kind="cubic",
            bounds_error=False,
            fill_value=(self.rho_l0_values[0], 0.0),
        )

    def _compute_rho_l0(self, a: float) -> float:
        """
        Compute the l-th multipole moment of the density at radius a.

        Parameters
        ----------
        a : float
            The spherical radius.

        Returns
        -------
        float
            The l-th multipole moment of the density at radius a.
        """

        def integrand(u: float) -> float:
            """Integrand for the l-th multipole moment."""
            R = float(a * np.sqrt(1 - u**2))
            z = a * u
            legendre_val = np.polynomial.legendre.Legendre.basis(self.multipole_order)(u)
            return float(self.model.density(R, z) * np.real(legendre_val))

        cl = float(np.sqrt(np.pi * (2 * self.multipole_order + 1)))
        result, _ = quad(integrand, -1, 1)
        return result * cl

    def __call__(self, a: float) -> float:
        """
        Evaluate the interpolated rho_l0 at radius a.

        Parameters
        ----------
        a : float
            The spherical radius.

        Returns
        -------
        float
            The interpolated l-th multipole moment of the density at radius a.
        """
        if a < self.a_min or a > self.a_max:
            return self._compute_rho_l0(a)
        return float(self.interpolator(a))


class MultipoleExpansion:
    """Class for multipole expansion of a galaxy model."""

    model: Model
    l_max: int

    def __init__(self, model: Model, l_max: int, /, *, Ngrid: int = 100) -> None:
        """
        Initialize the multipole expansion with a given model and maximum multipole order.

        Parameters
        ----------
        model : Model
            The galaxy density model.
        l_max : int
            Maximum multipole order.
        Ngrid : int, optional
            Number of grid points for interpolators, by default 100.
        """
        self.model = model
        self.l_max = l_max

        # Initialize interpolators for each multipole order
        self.rho_l0_interpolators = {}
        for multipole_order in range(0, l_max + 1, 2):
            self.rho_l0_interpolators[multipole_order] = RhoL0Interpolator(
                model,
                multipole_order,
                Ngrid=Ngrid,
            )
            logger.info("RhoL0Interpolator initialized for multipole order %d", multipole_order)

    def compute_potential(self, R: float, z: float) -> float:
        """
        Compute the gravitational potential at cylindrical coordinates (R, z) using the multipole expansion.

        Parameters
        ----------
        R : float
            The cylindrical radius.
        z : float
            The height above the plane.

        Returns
        -------
        float
            The gravitational potential at the given coordinates.
        """
        r = np.sqrt(R**2 + z**2)
        theta = np.arctan2(R, z)

        def external_integrand(a: float, multipole_order: int) -> float:
            """Auxiliary integrand for the external region."""
            return self.rho_l0_interpolators[multipole_order](a) * a ** (multipole_order + 2)

        def internal_integrand(a: float, multipole_order: int) -> float:
            """Auxiliary integrand for the internal region."""
            return self.rho_l0_interpolators[multipole_order](a) / a ** (multipole_order - 1)

        phi = 0.0
        for l in range(0, self.l_max + 1, 2):  # noqa: E741
            cl = np.sqrt((2 * l + 1) / (4 * np.pi))
            Pl = np.polynomial.legendre.Legendre.basis(l)(np.cos(theta))

            phi += (cl * Pl / (2 * l + 1)) * (
                quad(external_integrand, 0, r, args=(l,))[0] / r ** (l + 1)
                + quad(internal_integrand, r, np.inf, args=(l,))[0] * r**l
            )

        return float(-4 * np.pi * phi)


@register_asset("multipole_expansion_satoh.png")
def plot_multipole_expansion_satoh(path: Path | None = None) -> None:
    """Plot an example of the multipole expansion against the exact potential."""
    fig: Figure
    axs: Axes
    fig, axs = plt.subplots(
        1,
        1,
        figsize=(5, 5),
        sharex=True,
    )

    fig.subplots_adjust(
        left=0.1,
        right=0.92,
        bottom=0.1,
        top=0.92,
        wspace=0.05,
        hspace=0.05,
    )

    axs.set_xlabel(r"$R/a$")
    axs.set_ylabel(r"$z/a$")
    axs.set_xlim(0, 5)
    axs.set_ylim(0, 5)

    # configure ticks
    axs.xaxis.set_major_locator(MultipleLocator(1))
    axs.yaxis.set_major_locator(MultipleLocator(1))

    axs.xaxis.set_minor_locator(AutoMinorLocator())
    axs.yaxis.set_minor_locator(AutoMinorLocator())

    axs.tick_params(which="minor", length=3, color="gray", direction="in")
    axs.tick_params(which="major", length=6, direction="in")
    axs.tick_params(top=True, right=True, which="both")

    mdl = SatohModel(SatohModelParams(q=0.6))

    r = np.linspace(0.01, 5.0, 30)
    z = np.linspace(0.01, 5.0, 30)
    R, Z = np.meshgrid(r, z)
    levels = -1 / np.arange(1, 6, 0.7)

    # contours of exact potential
    P_exact = np.array([[mdl.potential(ri, zi) for ri in r] for zi in z])
    _ = axs.contour(
        R,
        Z,
        P_exact,
        levels=levels,
        colors="black",
        linestyles="solid",
        linewidths=1.5,
    )

    # contours of multipole expansion potential
    expansion = MultipoleExpansion(mdl, 4, Ngrid=500)
    P_expansion = np.array([[expansion.compute_potential(ri, zi) for ri in r] for zi in z])
    _ = axs.contour(
        R,
        Z,
        P_expansion,
        levels=levels,
        colors="black",
        linestyles="dashed",
        linewidths=1.5,
    )

    expansion = MultipoleExpansion(mdl, 6, Ngrid=500)
    P_expansion = np.array([[expansion.compute_potential(ri, zi) for ri in r] for zi in z])
    _ = axs.contour(
        R,
        Z,
        P_expansion,
        levels=levels,
        colors="black",
        linestyles="dotted",
        linewidths=1.5,
    )

    # Create legend using proxy artists
    legend_elements = [
        Line2D([0], [0], color="black", linestyle="solid", linewidth=1.5, label=r"Exact $b/a=0.6$"),
        Line2D([0], [0], color="black", linestyle="dashed", linewidth=1.5, label=r"$l_{\rm max} = 4$"),
        Line2D([0], [0], color="black", linestyle="dotted", linewidth=1.5, label=r"$l_{\rm max} = 6$"),
    ]
    axs.legend(handles=legend_elements, loc="upper right", frameon=False)

    if path:
        save_figure_if_changed(
            fig,
            path,
            dpi=150,
            bbox_inches="tight",
            facecolor="white",
            edgecolor="none",
        )
    else:
        plt.show()
