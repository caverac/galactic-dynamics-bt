"""Bifurcation analysis and orbit classification in logarithmic potentials.

This module provides tools for studying orbital structure in 2D logarithmic
potentials, including surface of section computation and loop/box orbit
classification.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
import logging
from pathlib import Path
from typing import Generic, TypeVar

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import numpy.typing as npt
from scipy.integrate import RK45


from galactic_dynamics_bt.chapter03.symplectic_integrators import SymplecticIntegrator
from galactic_dynamics_bt.utils.assets import register_asset, save_figure_if_changed


T = TypeVar("T")
logger = logging.getLogger(__name__)


class Model(ABC, Generic[T]):
    """Abstract base class for 2D potential models in Cartesian coordinates.

    Provides orbit integration, surface of section computation, and orbit
    classification for potentials defined in (x, y) coordinates.

    Parameters
    ----------
    params : T
        Model-specific parameters.
    """

    params: T

    def __init__(self, params: T) -> None:
        """Initialize the model with given parameters.

        Parameters
        ----------
        params : T
            Model-specific parameters.
        """
        self.params = params

    @abstractmethod
    def potential(self, x: float, y: float) -> float:
        """Calculate the gravitational potential at (x, y).

        Parameters
        ----------
        x : float
            Horizontal position.
        y : float
            Vertical position.

        Returns
        -------
        float
            Potential value Phi(x, y).
        """

    @abstractmethod
    def potential_gradient(self, x: float, y: float) -> tuple[float, float]:
        """Calculate the gradient of the potential at (x, y).

        Parameters
        ----------
        x : float
            Horizontal position.
        y : float
            Vertical position.

        Returns
        -------
        tuple[float, float]
            Gradient components (dPhi/dx, dPhi/dy).
        """

    def find_initial_conditions(self, E: float, x: float, py: float) -> tuple[float, float, float, float]:
        """Find initial conditions for orbit integration.

        Given energy E, position x at y=0, and momentum py, compute px from
        energy conservation.

        Parameters
        ----------
        E : float
            Total orbital energy.
        x : float
            Initial x position.
        py : float
            Initial y-momentum.

        Returns
        -------
        tuple[float, float, float, float]
            Initial conditions (x0, y0, px0, py0).

        Raises
        ------
        ValueError
            If energy is below potential or no real solution exists.
        """
        V = self.potential(x, 0.0)
        if E < V:
            raise ValueError("Energy is less than potential at given position.")
        y0 = 0.0
        px_squared = 2 * (E - V) - py**2
        if px_squared < 0:
            raise ValueError("No real solution for px with given energy and position.")
        px0 = np.sqrt(px_squared)
        return x, y0, px0, py

    def find_initial_conditions_on_section(self, E: float, x: float, px: float) -> tuple[float, float, float, float]:
        """Find initial conditions starting on the surface of section (y=0).

        Given (x, px) on the surface of section, compute py from energy
        conservation with py > 0 to ensure upward crossing.

        Parameters
        ----------
        E : float
            Total orbital energy.
        x : float
            Initial x position on section.
        px : float
            Initial x-momentum on section.

        Returns
        -------
        tuple[float, float, float, float]
            Initial conditions (x, 0, px, py).

        Raises
        ------
        ValueError
            If no real solution exists for py.
        """
        V = self.potential(x, 0.0)
        py_squared = 2 * (E - V) - px**2
        if py_squared < 0:
            raise ValueError("No real solution for py with given energy and position.")
        py0 = np.sqrt(py_squared)
        return x, 0.0, px, py0

    def integrate(
        self, E: float, t: npt.NDArray[np.float64], x: float, py: float, order: int = 4
    ) -> npt.NDArray[np.float64]:
        """Integrate orbit using symplectic integration.

        Parameters
        ----------
        E : float
            Total orbital energy.
        t : npt.NDArray[np.float64]
            Time array for integration.
        x : float
            Initial x position (at y=0).
        py : float
            Initial y-momentum.
        order : int, optional
            Integrator order (1, 2, 4, 6, 8). Default is 4.

        Returns
        -------
        npt.NDArray[np.float64]
            Solution array of shape (len(t), 4) with columns [x, y, px, py].
        """
        x0, y0, px0, py0 = self.find_initial_conditions(E, x, py)

        integrator = SymplecticIntegrator(
            lambda p: p,
            lambda q: np.array(self.potential_gradient(q[0], q[1])),
            np.array([x0, y0, px0, py0]),
            t,
            order=order,
        )

        return integrator.integrate()

    def is_loop_orbit(self, E: float, t: npt.NDArray[np.float64], x: float, py: float, order: int = 4) -> bool:
        """Classify orbit as loop or box based on angular momentum.

        Loop orbits circulate around the center (Lz keeps the same sign).
        Box orbits oscillate back and forth (Lz changes sign).

        Parameters
        ----------
        E : float
            Total orbital energy.
        t : npt.NDArray[np.float64]
            Time array for integration.
        x : float
            Initial x position (at y=0).
        py : float
            Initial y-momentum.
        order : int, optional
            Integrator order. Default is 4.

        Returns
        -------
        bool
            True if loop orbit, False if box orbit.
        """
        z = self.integrate(E, t, x, py, order=order)
        x_vals = z[:, 0]
        y_vals = z[:, 1]
        px_vals = z[:, 2]
        py_vals = z[:, 3]

        # Angular momentum Lz = x * py - y * px
        Lz = x_vals * py_vals - y_vals * px_vals

        # Loop orbit: Lz stays positive or stays negative (no sign changes)
        # Box orbit: Lz changes sign
        return bool(np.all(Lz > 0) or np.all(Lz < 0))

    def surface_of_section(
        self,
        E: float,
        t: npt.NDArray[np.float64],
        x: float,
        py: float,
        order: int = 4,
    ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        """Compute Poincare surface of section at y=0 with py>0.

        Parameters
        ----------
        E : float
            Total orbital energy.
        t : npt.NDArray[np.float64]
            Time array for integration.
        x : float
            Initial x position.
        py : float
            Initial y-momentum.
        order : int, optional
            Integrator order. Default is 4.

        Returns
        -------
        tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]
            Arrays (x_crossings, px_crossings) at y=0.
        """
        z = self.integrate(E, t, x, py, order=order)
        return self._extract_crossings(z)

    def surface_of_section_from_section(
        self,
        E: float,
        t: npt.NDArray[np.float64],
        x0: float,
        px0: float,
        order: int = 4,
    ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        """Compute surface of section starting from a point on the section.

        Parameters
        ----------
        E : float
            Total orbital energy.
        t : npt.NDArray[np.float64]
            Time array for integration.
        x0 : float
            Initial x position on section.
        px0 : float
            Initial x-momentum on section.
        order : int, optional
            Integrator order. Default is 4.

        Returns
        -------
        tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]
            Arrays (x_crossings, px_crossings) at y=0.
        """
        x, y, px, py = self.find_initial_conditions_on_section(E, x0, px0)
        integrator = SymplecticIntegrator(
            lambda p: p,
            lambda q: np.array(self.potential_gradient(q[0], q[1])),
            np.array([x, y, px, py]),
            t,
            order=order,
        )
        z = integrator.integrate()
        return self._extract_crossings(z)

    def _extract_crossings(self, z: npt.NDArray[np.float64]) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        """Extract y=0 crossings with py>0 from integrated orbit.

        Parameters
        ----------
        z : npt.NDArray[np.float64]
            Integrated orbit array of shape (N, 4) with columns [x, y, px, py].

        Returns
        -------
        tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]
            Arrays (x_crossings, px_crossings) at y=0 with py>0.
        """
        x_vals = z[:, 0]
        y_vals = z[:, 1]
        px_vals = z[:, 2]
        py_vals = z[:, 3]

        crossings_x = []
        crossings_px = []

        for i in range(1, len(y_vals)):
            if y_vals[i - 1] < 0 <= y_vals[i] and py_vals[i] > 0:
                x_cross, px_cross = self.find_exact_intersection(x_vals[i], px_vals[i], y_vals[i], py_vals[i])
                crossings_x.append(x_cross)
                crossings_px.append(px_cross)

        return np.array(crossings_x), np.array(crossings_px)

    def find_exact_intersection(self, x: float, px: float, y: float, py: float) -> tuple[float, float]:
        """Find exact y=0 crossing via backward integration.

        Integrates backward from (x, y, px, py) to find the exact point
        where the orbit crosses y=0.

        Parameters
        ----------
        x : float
            Position x near crossing.
        px : float
            Momentum px near crossing.
        y : float
            Position y (should be small positive).
        py : float
            Momentum py (should be positive).

        Returns
        -------
        tuple[float, float]
            Exact (x, px) at y=0.
        """

        def derivatives(y_curr: float, z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
            _t_curr, x_curr, px_curr, py_curr = z
            dpx_dt, dpy_dt = self.potential_gradient(x_curr, y_curr)
            dt_dy = 1 / py_curr
            dx_dy = px_curr / py_curr
            dpx_dy = -dpx_dt / py_curr
            dpy_dy = -dpy_dt / py_curr
            return np.array([dt_dy, dx_dy, dpx_dy, dpy_dy])

        z = np.array([0, x, px, py])
        sol = RK45(derivatives, y, z, t_bound=0, max_step=1e-3)
        while sol.status == "running":
            sol.step()
        _t_final, x_final, px_final, _py_final = sol.y

        return x_final, px_final


@dataclass
class LogarithmicPotentialParams:
    """Parameters for the logarithmic potential.

    Attributes
    ----------
    q : float
        Flattening parameter (0 < q <= 1). q=1 is spherical.
    v0 : float
        Characteristic velocity scale.
    Rc : float
        Core radius for softening.
    """

    q: float
    v0: float
    Rc: float


class LogarithmicPotential(Model[LogarithmicPotentialParams]):
    """Logarithmic potential model: Phi = (v0^2/2) * ln(x^2 + y^2/q^2 + Rc^2).

    This potential produces flat rotation curves and supports both loop
    and box orbits depending on the flattening parameter q.
    """

    def potential(self, x: float, y: float) -> float:
        """Calculate the logarithmic potential at (x, y).

        Parameters
        ----------
        x : float
            Horizontal position.
        y : float
            Vertical position.

        Returns
        -------
        float
            Potential value.
        """
        v0 = self.params.v0
        q = self.params.q
        Rc = self.params.Rc
        return 0.5 * v0**2 * float(np.log(x**2 + (y**2) / (q**2) + Rc**2))

    def potential_gradient(self, x: float, y: float) -> tuple[float, float]:
        """Calculate the gradient of the logarithmic potential.

        Parameters
        ----------
        x : float
            Horizontal position.
        y : float
            Vertical position.

        Returns
        -------
        tuple[float, float]
            Gradient components (dPhi/dx, dPhi/dy).
        """
        q = self.params.q
        v0 = self.params.v0
        Rc = self.params.Rc
        denom = x**2 + (y**2) / (q**2) + Rc**2
        dR = (v0**2) * x / denom
        dz = (v0**2) * y / (q**2 * denom)
        return dR, dz


class AxesArray:
    """
    Helper class for type hinting 2D arrays of matplotlib Axes.

    Provides proper type hints for subplot arrays returned by
    plt.subplots() with multiple rows and columns. This ensures
    type checkers understand that axs[i, j] returns an Axes object.
    """

    def __getitem__(self, _index: tuple[int, int] | int) -> Axes:  # type: ignore[empty-body]
        """Return Axes object at given 2D index."""
        ...


@register_asset("surfaces_of_section.png")
def plot_surfaces_of_section(path: Path | None = None) -> None:
    """Create 2x2 plot of surfaces of section for different q values.

    Generates a figure with four panels showing Poincare surfaces of section
    for q = 0.7, 0.8, 0.9, 1.0. Each panel shows orbits sampled at different
    initial conditions along with the zero velocity curve boundary.

    Parameters
    ----------
    path : Path | None, optional
        File path to save the figure. If None, displays interactively.
    """
    Rc = 0.14
    t: npt.NDArray[np.float64] = np.arange(0, 200, 0.01, dtype=np.float64)
    q_values = [0.7, 0.8, 0.9, 1.0]

    fig: Figure
    axs: AxesArray

    fig, axs = plt.subplots(2, 2, figsize=(6, 5), sharex=True, sharey=True)
    fig.subplots_adjust(
        left=0.12,
        right=0.98,
        bottom=0.1,
        top=0.92,
        wspace=0.08,
        hspace=0.08,
    )

    for idx, q in enumerate(q_values):

        row, col = idx // 2, idx % 2
        ax: Axes = axs[row, col]

        model = LogarithmicPotential(LogarithmicPotentialParams(q=q, v0=1.0, Rc=Rc))
        E = model.potential(5 * Rc, 0.0)

        # Sample initial conditions (x0, px0) on the section plane
        # Sample along px=0 at different x values (spans from center to edge)
        for x0 in np.linspace(0.05, 0.85, 8):
            try:
                x_sec, px_sec = model.surface_of_section_from_section(E, t, x0, 0.0, order=2)
                ax.plot(x_sec, px_sec, "k.", markersize=1)
                ax.plot(-x_sec, -px_sec, "k.", markersize=1)
            except ValueError:
                pass

        # Also sample at different px values along x=0.5 (vertical slice)
        for px0 in np.linspace(-0.8, 0.8, 8):
            try:
                x_sec, px_sec = model.surface_of_section_from_section(E, t, 0.5, px0, order=2)
                ax.plot(x_sec, px_sec, "k.", markersize=1)
                ax.plot(-x_sec, -px_sec, "k.", markersize=1)
            except ValueError:
                pass

        # Zero velocity curve: at y=0, E = V(x,0) + 0.5*px^2
        # So px = +/- sqrt(2*(E - V(x,0)))
        x_zvc = np.linspace(-1.0, 1.0, 200)
        V_zvc = np.array([model.potential(xi, 0.0) for xi in x_zvc])
        px_squared = 2 * (E - V_zvc)
        valid = px_squared >= 0
        ax.plot(x_zvc[valid], np.sqrt(px_squared[valid]), "k-", linewidth=1)
        ax.plot(x_zvc[valid], -np.sqrt(px_squared[valid]), "k-", linewidth=1)

        # Formatting
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(which="minor", length=3, color="gray", direction="in")
        ax.tick_params(which="major", length=6, direction="in")
        ax.tick_params(top=True, right=True, which="both")
        ax.set_xlim(-1.0, 1.0)
        ax.set_ylim(-2.0, 2.0)
        ax.text(0.95, 0.95, f"$q={q}$", transform=ax.transAxes, ha="right", va="top", fontsize=10)

        if col == 0:
            ax.set_ylabel(r"$p_x$")
        if row == 1:
            ax.set_xlabel(r"$x$")

    if path:
        save_figure_if_changed(fig, path, dpi=150, bbox_inches="tight", facecolor="white", edgecolor="none")
    else:
        plt.show()


@register_asset("loop_orbit_fraction.png")
def plot_loop_orbit_fraction(path: Path | None = None) -> None:
    """Plot fraction of loop orbits as a function of flattening parameter q.

    Samples a grid of initial conditions and classifies each orbit as loop
    or box based on angular momentum conservation. Shows how the orbital
    structure transitions from box-dominated (low q) to loop-dominated (q=1).

    Parameters
    ----------
    path : Path | None, optional
        File path to save the figure. If None, displays interactively.
    """
    Rc = 0.14
    t: npt.NDArray[np.float64] = np.arange(0, 100, 0.05, dtype=np.float64)

    q_values = [0.7, 0.8, 0.85, 0.9, 0.95, 0.99, 1.0]
    x0_values = np.linspace(0.1, 0.5, 12)
    py0_values = np.linspace(0.4, 1.0, 12)

    loop_fractions = []

    for q in q_values:
        model = LogarithmicPotential(LogarithmicPotentialParams(q=q, v0=1.0, Rc=Rc))
        E = model.potential(10 * Rc, 0.0)

        loop_count = 0
        total_count = 0

        for x0 in x0_values:
            for py0 in py0_values:
                try:
                    is_loop = model.is_loop_orbit(E, t, x0, py0, order=2)
                    if is_loop:
                        loop_count += 1
                    total_count += 1
                except ValueError:
                    pass

        fraction = loop_count / total_count if total_count > 0 else 0.0
        loop_fractions.append(fraction)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(q_values, loop_fractions, "o-", markersize=8, linewidth=2, color="k")

    ax.set_xlabel(r"$q$")
    ax.set_ylabel("Fraction of loop orbits")
    ax.set_xlim(0.65, 1.05)
    ax.set_ylim(-0.05, 1.05)

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="minor", length=3, color="gray", direction="in")
    ax.tick_params(which="major", length=6, direction="in")
    ax.tick_params(top=True, right=True, which="both")

    if path:
        save_figure_if_changed(fig, path, dpi=150, bbox_inches="tight", facecolor="white", edgecolor="none")
    else:
        plt.show()


if __name__ == "__main__":  # pragma: no cover
    plot_loop_orbit_fraction()
    # plot_surfaces_of_section()
