"""Surface of section analysis for stellar orbits in logarithmic potentials.

This module provides some functionality for studying stellar orbits in galaxies
with logarithmic potentials. It includes:

- Abstract base class for galactic potential models
- Logarithmic potential implementation with flattening parameter
- Symplectic orbit integration for long-term stability
- Surface of section analysis for studying orbital structure
- Visualization tools for orbits and Poincare sections

The logarithmic potential is commonly used in galactic dynamics as it provides
a simple yet realistic model for galaxy halos and rotation curves. The form is:

    Phi(R,z) = (v0^2/2) * ln(R^2 + z^2/q^2)

where v0 sets the velocity scale and q is the flattening parameter.

Examples
--------
Basic orbit integration:

>>> params = LogarithmicPotentialParams(q=0.9, v0=1.0, Rc=0.0)
>>> model = LogarithmicPotential(params, Lz=0.2)
>>> t, R, z, pR, pz = integrate_orbit(model, E=-0.8, R0=0.35, pR0=0.1, t_max=5.0, dt=0.01)

Surface of section analysis:

>>> R_sec, pR_sec = surface_of_section(model, R, z, pR, pz)

"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
import logging
from pathlib import Path
from typing import Generic, TypeVar

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import numpy as np
from scipy.integrate import RK45
from scipy.optimize import root_scalar

from galactic_dynamics_bt.chapter03.symplectic_integrators import SymplecticIntegrator


T = TypeVar("T")
logger = logging.getLogger(__name__)


class Model(ABC, Generic[T]):
    """
    Abstract base class for galaxy density and potential models.

    Provides the interface for galactic potential models used in orbit integration.
    Implements common functionality like effective potential calculation that
    includes centrifugal force terms due to angular momentum conservation.

    The class uses cylindrical coordinates (R, z, phi) where R is the radial
    distance in the galactic plane, z is the height above/below the plane,
    and phi is the azimuthal angle. Angular momentum Lz is conserved.

    Parameters
    ----------
    params : T
        Model-specific parameters (e.g., mass scales, length scales)
    Lz : float
        Specific angular momentum in the z-direction (conserved quantity)

    Attributes
    ----------
    params : T
        The model parameters
    Lz : float
        Specific angular momentum
    """

    params: T
    Lz: float

    def __init__(self, params: T, Lz: float) -> None:
        """
        Initialize the model with given parameters.

        Parameters
        ----------
        params : T
            Model-specific parameters object containing physical constants
            and scaling factors needed to define the potential.
        Lz : float
            Specific angular momentum about the z-axis. This quantity is
            conserved in axisymmetric potentials and determines the
            centrifugal barrier in the effective potential.

        Examples
        --------
        >>> params = LogarithmicPotentialParams(q=0.9, v0=1.0, Rc=0.0)
        >>> model = LogarithmicPotential(params, Lz=0.2)
        """
        self.params = params
        self.Lz = Lz

    @abstractmethod
    def potential(self, R: float, z: float) -> float:
        """
        Calculate the gravitational potential at cylindrical coordinates (R, z).

        Parameters
        ----------
        R : float
            Radial distance from the z-axis (R ≥ 0)
        z : float
            Height above/below the galactic plane

        Returns
        -------
        float
            Gravitational potential Φ(R,z) in model units

        Notes
        -----
        The potential should be defined such that the force is F = -grad Phi.
        For bound orbits, the potential is typically negative with
        Phi → 0 as r → ∞.
        """

    @abstractmethod
    def potential_gradient(self, R: float, z: float) -> tuple[float, float]:
        """
        Calculate the R, z components of the gradient of the gravitational potential.

        Parameters
        ----------
        R : float
            Radial distance from the z-axis (R > 0)
        z : float
            Height above/below the galactic plane

        Returns
        -------
        tuple[float, float]
            Gradient components (dPhi/dR, dPhi/dz)

        Notes
        -----
        Returns the gradient of the potential.
        """

    @abstractmethod
    def zero_velocity_curve(self, E: float, Lz: float) -> tuple[np.ndarray, np.ndarray]:
        """
        Calculate the zero-velocity curve for given energy and angular momentum.

        The zero-velocity curve defines the boundary in configuration space
        where the kinetic energy becomes zero: E = Phi_eff(R,z).

        Parameters
        ----------
        E : float
            Total orbital energy
        Lz : float
            Angular momentum about z-axis

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            Arrays (R, z) defining the zero-velocity curve boundary

        Notes
        -----
        The particle cannot access regions where E < Phi_eff(R,z).
        This curve is essential for understanding the allowed
        orbital regions and for surface of section analysis.
        """

    def effective_potential(self, R: float, z: float) -> float:
        """
        Calculate the effective potential including centrifugal term.

        The effective potential combines the gravitational potential with
        the centrifugal potential due to angular momentum conservation:

        Phi_eff(R,z) = Phi(R,z) + Lz^2/(2R^2)

        Parameters
        ----------
        R : float
            Radial distance from z-axis (R > 0)
        z : float
            Height above/below galactic plane

        Returns
        -------
        float
            Effective potential energy

        Notes
        -----
        The centrifugal term Lz^2/(2R^2) creates a repulsive barrier
        that prevents particles from reaching R = 0. This reduces
        the 3D orbital problem to an effective 2D problem in the
        meridional plane (R, z).
        """
        return self.potential(R, z) + 0.5 * (self.Lz**2) / (R**2)

    def effective_potential_gradient(self, R: float, z: float) -> tuple[float, float]:
        """
        Calculate the gradient of the effective potential including centrifugal term.

        Computes the negative gradient of the effective potential:
        dPhi_eff/dR = dPhi/dR - Lz^2/R^3
        dPhi_eff/dz = dPhi/dz

        Parameters
        ----------
        R : float
            Radial distance from z-axis (R > 0)
        z : float
            Height above/below galactic plane

        Returns
        -------
        tuple[float, float]
            Gradient components (dPhi_eff/dR, dPhi_eff/dz)

        Notes
        -----
        The radial component includes the centrifugal force Lz^2/R^3
        which acts outward, opposing gravitational attraction.
        The vertical component is unchanged from the gravitational force.
        """
        dR, dZ = self.potential_gradient(R, z)
        dR -= self.Lz**2 / R**3
        return dR, dZ

    def find_initial_conditions(self, E: float, R: float, pR: float) -> tuple[float, float, float, float]:
        """
        Find initial conditions for orbit integration given energy and starting point.

        For a given total energy E and starting position (R, z=0) with radial
        momentum pR, calculates the required vertical momentum pz to satisfy
        the energy constraint: E = T + Phi_eff where T = (pR^2 + pz^2)/2.

        Parameters
        ----------
        E : float
            Total orbital energy
        R : float
            Initial radial position (R > 0)
        pR : float
            Initial radial momentum

        Returns
        -------
        tuple[float, float, float, float]
            Initial conditions (R, z, pR, pz) where z=0 and pz is determined
            from energy conservation

        Notes
        -----
        Assumes the orbit starts at the galactic plane (z = 0) and calculates
        the required vertical momentum from energy conservation:
        pz = √[2(E - pR^2/2 - Phi_eff(R,0))]
        """
        E -= 0.5 * pR**2
        z = 0.0
        pz = np.sqrt(2 * (E - self.effective_potential(R, z)))
        return R, z, pR, pz


@dataclass
class LogarithmicPotentialParams:
    """
    Parameters for the logarithmic potential model.

    The logarithmic potential is defined as:
    Phi(R,z) = (v0^2/2) * ln(R^2 + z^2/q^2)
    This potential produces flat rotation curves and is commonly used
    to model galaxy halos and dark matter distributions.

    Attributes
    ----------
    q : float
        Flattening parameter (0 < q ≤ 1). Values q < 1 correspond to
        oblate (flattened) systems, while q = 1 gives a spherical potential.
        Typical values for galaxy halos are q ~ 0.6-0.9.
    v0 : float
        Characteristic velocity scale. Sets the amplitude of the rotation
        curve and the overall velocity scale of the system. In physical
        units, this would typically be in km/s.
    Rc : float
        Core radius parameter. Introduces a softening length scale that
        prevents the potential from diverging at the center (R=0, z=0).
        Useful for modeling systems with a finite central density.

    Examples
    --------
    >>> # Nearly spherical halo
    >>> params = LogarithmicPotentialParams(q=0.9, v0=220.0, Rc=0.0)  # 220 km/s
    >>>
    >>> # Significantly flattened system
    >>> params = LogarithmicPotentialParams(q=0.6, v0=1.0, Rc=0.0)   # Normalized units
    """

    q: float
    v0: float
    Rc: float


class LogarithmicPotential(Model[LogarithmicPotentialParams]):
    """
    Logarithmic potential model for galactic dynamics.

    Implements the logarithmic potential commonly used in galactic dynamics:
    Phi(R,z) = (v0^2/2) * ln(R^2 + z^2/q^2)

    This potential has several attractive properties:
    - Produces flat rotation curves: v_c(∞) = v0
    - Simple analytical form allowing exact calculations
    - Realistic for modeling galaxy halos and dark matter
    - Supports both bound and unbound orbits depending on energy

    The flattening parameter q allows modeling of oblate systems
    typical of galaxy halos formed through hierarchical structure formation.

    Parameters
    ----------
    params : LogarithmicPotentialParams
        Container with potential parameters (q, v0)
    Lz : float
        Specific angular momentum about z-axis

    Examples
    --------
    Create a model and compute properties:

    >>> params = LogarithmicPotentialParams(q=0.8, v0=1.0, Rc=0.0)
    >>> model = LogarithmicPotential(params, Lz=0.2)
    >>> phi = model.potential(1.0, 0.5)  # Potential at (R,z) = (1,0.5)
    >>> F_R, F_z = model.potential_gradient(1.0, 0.5)  # Forces

    Compute circular velocity at radius R:

    >>> R = 1.0
    >>> dPhi_dR, _ = model.potential_gradient(R, 0.0)
    >>> v_circ = np.sqrt(R * dPhi_dR)  # Should equal v0 for large R
    """

    def potential(self, R: float, z: float) -> float:
        """
        Calculate the logarithmic potential at (R, z).

        Implements: Phi(R,z) = (v0^2/2) * ln(R^2 + z^2/q^2 + Rc^2)

        Parameters
        ----------
        R : float
            Radial coordinate (R > 0)
        z : float
            Vertical coordinate

        Returns
        -------
        float
            Gravitational potential value

        Examples
        --------
        >>> model = LogarithmicPotential(LogarithmicPotentialParams(q=0.9, v0=1.0, Rc=0.0), Lz=0.2)
        >>> phi_center = model.potential(1.0, 0.0)    # At galactic plane
        >>> phi_above = model.potential(1.0, 0.5)     # Above plane
        >>> print(f"Potential difference: {phi_above - phi_center:.3f}")
        Potential difference: 0.134
        """
        q = self.params.q
        v0 = self.params.v0
        Rc = self.params.Rc
        return 0.5 * v0**2 * float(np.log(R**2 + (z**2) / (q**2) + Rc**2))

    def potential_gradient(self, R: float, z: float) -> tuple[float, float]:
        """
        Calculate the gradient for logarithmic potential.

        Computes the negative gradient: F = -grad Phi

        dPhi/dR = -v0^2 * R / (R^2 + z^2/q^2 + Rc^2)
        dPhi/dz = -v0^2 * z / (q^2(R^2 + z^2/q^2 + Rc^2))

        Parameters
        ----------
        R : float
            Radial coordinate (R > 0)
        z : float
            Vertical coordinate

        Returns
        -------
        tuple[float, float]
            Gradient components (-dPhi/dR, -dPhi/dz)
        """
        q = self.params.q
        v0 = self.params.v0
        Rc = self.params.Rc
        denom = R**2 + (z**2) / (q**2) + Rc**2
        dR = (v0**2) * R / denom
        dz = (v0**2) * z / (q**2 * denom)
        return dR, dz

    def zero_velocity_curve(self, E: float, Lz: float) -> tuple[np.ndarray, np.ndarray]:
        """
        Calculate the zero-velocity curve for given energy and angular momentum.

        Finds the boundary where kinetic energy vanishes: E = Phi_eff(R,z).
        For the logarithmic potential, this defines the allowed orbital region.

        Parameters
        ----------
        E : float
            Total orbital energy
        Lz : float
            Angular momentum about z-axis

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            Arrays (R, z) defining the zero-velocity curve

        Notes
        -----
        The zero-velocity curve separates allowed (E > Phi_eff) from
        forbidden (E < Phi_eff) regions. At this boundary:

        E = Phi(R,z) + Lz^2/(2R^2)

        For bound orbits (E < 0), this creates a finite allowed region.
        For unbound orbits (E ≥ 0), particles can escape to infinity.

        The calculation first finds the radial turning points where
        the effective potential in the plane equals the energy, then
        solves for z(R) at each radius.

        """

        def f(x: float) -> float:
            return float(E - self.params.v0**2 * 0.5 * np.log(x**2) - 0.5 * (Lz**2) / (x**2))

        Rmin = root_scalar(f, bracket=(0.01, 0.3), method="brentq").root
        Rmax = root_scalar(f, bracket=(0.3, 2.0), method="brentq").root

        R = np.linspace(Rmin, Rmax, 500, endpoint=True)
        return R, np.sqrt(
            2
            * (E - 0.5 * self.params.v0**2 * np.log(R**2 + self.params.Rc**2) - 0.5 * (Lz**2) / (R**2))
            * (self.params.q**2)
        )


def integrate_orbit(
    model: Model[T],
    E: float,
    R0: float,
    pR0: float,
    t_max: float,
    dt: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Integrate an orbit in the given model.

    Parameters
    ----------
    model : Model[T]
        The galaxy model to use for the orbit integration.
    E : float
        The total energy of the orbit.
    R0 : float
        The initial radial coordinate.
    pR0 : float
        The initial radial momentum.
    t_max : float
        The maximum integration time.
    dt : float
        The time step for integration.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        Arrays of time, R, z, pR, and pz values along the orbit.
    """
    R0, z0, pR0, pz0 = model.find_initial_conditions(E, R0, pR0)

    t = np.arange(0, t_max, dt)
    integrator = SymplecticIntegrator(
        lambda p: p,
        lambda q: np.array(model.effective_potential_gradient(q[0], q[1])),
        np.array([R0, z0, pR0, pz0]),
        t,
        order=2,
        H=lambda q, p: model.effective_potential(q[0], q[1]) + 0.5 * (p[0] ** 2 + p[1] ** 2),
    )
    y = integrator.integrate()
    R = y[:, 0]
    z = y[:, 1]
    pR = y[:, 2]
    pz = y[:, 3]

    return t, R, z, pR, pz


def find_exact_intersection(R: float, z: float, pR: float, pz: float, model: Model[T]) -> tuple[float, float]:
    """
    Find exact intersection with Poincaré section plane using numerical integration.

    Uses backward integration from a point (R, z, pR, pz) to find the exact
    location where the orbit crosses z = 0 with pz > 0. This provides accurate
    surface of section coordinates by integrating the equations of motion
    backward to the crossing plane.

    Parameters
    ----------
    R : float
        Current radial coordinate
    z : float
        Current vertical coordinate (should be > 0)
    pR : float
        Current radial momentum
    pz : float
        Current vertical momentum (should be > 0)
    model : Model[T]
        Galaxy potential model for force calculations

    Returns
    -------
    tuple[float, float]
        Exact coordinates (R, pR) at the z = 0 crossing

    Notes
    -----
    The integration uses z as the independent variable and integrates
    backward from z > 0 to z = 0. The equations of motion in this
    parameterization are:

    dt/dz = 1/pz
    dpz/dz = -(dPhi/dz)/pz
    dR/dz = pR/pz
    dpR/dz = -(dPhi/dR)/pz

    This method is essential for accurate surface of section analysis
    since linear interpolation between grid points introduces errors
    in the Poincare coordinates.
    """

    def derivatives(z_curr: float, y: np.ndarray) -> np.ndarray:
        t_curr, pz_curr, R_curr, pR_curr = y  # pylint: disable=unused-variable
        dR_dt, dz_dt = model.effective_potential_gradient(R_curr, z_curr)
        dpR_dz = -dR_dt / pz_curr
        dpz_dz = -dz_dt / pz_curr
        dR_dz = pR_curr / pz_curr
        return np.array([1 / pz_curr, dpz_dz, dR_dz, dpR_dz])

    y = np.array([0, pz, R, pR])
    sol = RK45(derivatives, z, y, t_bound=0, max_step=1e-3)
    while sol.status == "running":
        sol.step()
    t_final, pz_final, R_final, pR_final = sol.y  # pylint: disable=unused-variable

    return R_final, pR_final


def surface_of_section(
    model: Model[T],
    R: np.ndarray,
    z: np.ndarray,
    pR: np.ndarray,
    pz: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute Poincaré surface of section for orbital analysis.

    Generates a surface of section by recording (R, pR) coordinates
    whenever the orbit crosses the z = 0 plane with pz > 0. This
    technique reduces the 4D phase space to a 2D map, revealing
    the underlying structure of orbital motion.

    Parameters
    ----------
    model : Model[T]
        Galaxy potential model for exact intersection calculations
    R : np.ndarray
        Array of radial coordinates from orbit integration
    z : np.ndarray
        Array of vertical coordinates from orbit integration
    pR : np.ndarray
        Array of radial momenta from orbit integration
    pz : np.ndarray
        Array of vertical momenta from orbit integration

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Arrays of (R, pR) coordinates at each z = 0 crossing with pz > 0

    Notes
    -----
    The surface of section is constructed by monitoring sign changes
    in the z coordinate and recording crossings where:
    - z changes from negative to positive (z[i-1] < 0 <= z[i])
    - The vertical momentum is positive (pz > 0)

    Each crossing is refined using exact numerical integration to
    eliminate discretization errors from the orbit integration grid.
    This ensures accurate representation of invariant curves and
    chaotic regions in phase space.

    For integrable systems, orbits appear as closed curves or
    discrete points on the surface of section. Chaotic orbits
    fill regions of phase space with apparently random points.

    """
    crossings_R = []
    crossings_pR = []

    for i in range(1, len(z)):
        if z[i - 1] < 0 <= z[i] and pz[i] > 0:

            R_exact, pR_exact = find_exact_intersection(R[i], z[i], pR[i], pz[i], model)
            crossings_R.append(R_exact)
            crossings_pR.append(pR_exact)

    return np.array(crossings_R), np.array(crossings_pR)


class AxesArray:
    """
    Helper class for type hinting 2D arrays of matplotlib Axes.

    Provides proper type hints for subplot arrays returned by
    plt.subplots() with multiple rows and columns. This ensures
    type checkers understand that axs[i, j] returns an Axes object.
    """

    def __getitem__(self, index: tuple[int, int] | int) -> Axes:  # type: ignore[empty-body]
        """Return Axes object at given 2D index."""
        ...


def plot_orbit(q: float, axs: Axes) -> None:
    """
    Plot a single orbit in the meridional (R, z) plane.

    Integrates and visualizes one orbital trajectory in the meridional
    plane for a logarithmic potential with specified flattening. The
    orbit is computed using symplectic integration for long-term stability.

    Parameters
    ----------
    q : float
        Flattening parameter for the logarithmic potential (0 < q ≤ 1)
    axs : Axes
        Matplotlib axes object for plotting

    Notes
    -----
    Fixed parameters used:
    - v0 = 1.0 (velocity scale)
    - Lz = 0.2 (angular momentum)
    - E = -0.8 (total energy)
    - Initial conditions: R0=0.35, pR0=0.1
    - Integration: t_max=20.0, dt=0.01

    The orbit shows the projection of 3D motion onto the meridional
    plane (R, z). Due to angular momentum conservation, the motion
    is confined to this 2D plane in the effective potential.
    """
    model = LogarithmicPotential(LogarithmicPotentialParams(q=q, v0=1.0, Rc=0.0), Lz=0.2)
    t, R, z, pR, pz = integrate_orbit(  # pylint: disable=unused-variable
        model,
        E=-0.8,
        R0=0.35,
        pR0=0.1,
        t_max=20.0,
        dt=0.01,
    )
    axs.xaxis.set_minor_locator(AutoMinorLocator())
    axs.yaxis.set_minor_locator(AutoMinorLocator())

    axs.set_xlim(0, 0.5)
    axs.set_ylim(-0.2, 0.2)

    axs.tick_params(which="minor", length=3, color="gray", direction="in")
    axs.tick_params(which="major", length=6, direction="in")
    axs.tick_params(top=True, right=True, which="both")

    axs.plot(R, z, "k", linewidth=1)


# pylint: disable=too-many-locals
def plot_surface_of_section(q: float, axs: Axes) -> None:
    """
    Plot Poincaré surface of section with additional reference curves.

    Creates a comprehensive surface of section plot showing:
    1. Orbit points at z=0 crossings with pz>0
    2. Zero-velocity curves (boundaries of allowed motion)
    3. Constant total angular momentum curves (if Lz were conserved)

    Parameters
    ----------
    q : float
        Flattening parameter for the logarithmic potential (0 < q ≤ 1)
    axs : Axes
        Matplotlib axes object for plotting

    Notes
    -----
    Fixed parameters used:
    - v0 = 1.0 (velocity scale)
    - Lz = 0.2 (specific angular momentum)
    - E = -0.8 (total energy)
    - Initial conditions: R0=0.35, pR0=0.1
    - Integration: t_max=100.0, dt=0.01 (longer than orbit plot)

    Reference curves:
    - Solid black lines: Zero-velocity curves at z = ±z_max(R)
    - Dashed black lines: Constant total angular momentum curves

    The total angular momentum L includes both the z-component Lz
    and contributions from meridional motion:
    L = √[(z*pR - R*pz)^2 + (R^2 + z^2)*Lz^2/R^2]
    """
    model = LogarithmicPotential(LogarithmicPotentialParams(q=q, v0=1.0, Rc=0.0), Lz=0.2)
    t, R, z, pR, pz = integrate_orbit(  # pylint: disable=unused-variable
        model,
        E=-0.8,
        R0=0.35,
        pR0=0.1,
        t_max=100.0,
        dt=0.01,
    )
    axs.yaxis.set_major_locator(MultipleLocator(0.5))

    axs.xaxis.set_minor_locator(AutoMinorLocator())
    axs.yaxis.set_minor_locator(AutoMinorLocator())

    axs.set_xlim(0, 0.5)
    axs.set_ylim(-0.8, 0.8)

    axs.tick_params(which="minor", length=3, color="gray", direction="in")
    axs.tick_params(which="major", length=6, direction="in")
    axs.tick_params(top=True, right=True, which="both")

    R_sec, pR_sec = surface_of_section(model, R, z, pR, pz)
    axs.plot(R_sec, pR_sec, "k.", markersize=1)

    # zero-velocity curve
    E = -0.8
    R_zvc, z_zvc = model.zero_velocity_curve(E, model.Lz)
    axs.plot(R_zvc, z_zvc, "k", linewidth=1)
    axs.plot(R_zvc, -z_zvc, "k", linewidth=1)

    # conserved total angular momentum curve
    L = np.sqrt((z[0] * pR[0] - R[0] * pz[0]) ** 2 + (R[0] ** 2 + z[0] ** 2) * model.Lz**2 / R[0] ** 2)

    R_tac, pR_tac = model.zero_velocity_curve(E, float(L))
    axs.plot(R_tac, pR_tac, "k--", linewidth=1)
    axs.plot(R_tac, -pR_tac, "k--", linewidth=1)


def plot_orbital_properties(path: Path | None = None) -> None:
    """
    Create comprehensive 2x2 plot comparing orbital properties for different flattenings.

    Generates a multi-panel figure showing the effect of potential flattening
    on orbital structure. Compares two flattening parameters (q=0.9 and q=0.6)
    through both orbital trajectories and surface of section analysis.

    Panel Layout:
    - Top row: Orbital trajectories in (R, z) plane
      - Left (q=0.9): Nearly spherical potential
      - Right (q=0.6): Significantly flattened potential
    - Bottom row: Poincaré surfaces of section at z=0
      - Left (q=0.9): Surface of section for spherical case
      - Right (q=0.6): Surface of section for flattened case

    Returns
    -------
    None
        Displays the plot using plt.show()

    Notes
    -----
    The comparison reveals how potential flattening affects:
    - Orbital shapes and allowed regions
    - Phase space structure and integrability
    - Conservation properties and resonances

    More flattened potentials (smaller q) typically show:
    - Increased vertical confinement
    - More complex phase space structure
    - Potential for chaotic motion

    Figure formatting:
    - Figure size: 6x5 inches
    - Shared x-axis for vertical comparison
    - Optimized spacing for publication quality

    """
    fig: Figure
    axs: AxesArray
    fig, axs = plt.subplots(
        2,
        2,
        figsize=(6, 5),
        sharex=True,
    )

    fig.subplots_adjust(
        left=0.1,
        right=0.98,
        bottom=0.1,
        top=0.92,
        wspace=0.15,
        hspace=0.05,
    )

    # top left: orbit in (R, z)
    axs[0, 0].set_ylabel(r"$z$")
    plot_orbit(0.9, axs[0, 0])

    # top right: orbit in (R, z)
    plot_orbit(0.6, axs[0, 1])

    # bottom left: surface of section
    axs[1, 0].set_xlabel(r"$R$")
    axs[1, 0].set_ylabel(r"$p_R$")
    plot_surface_of_section(0.9, axs[1, 0])

    # bottom right: surface of section
    axs[1, 1].set_xlabel(r"$R$")
    plot_surface_of_section(0.6, axs[1, 1])

    if path:
        fig.savefig(
            path,
            dpi=150,
            bbox_inches="tight",
            facecolor="white",
            edgecolor="none",
        )
    else:
        plt.show()


def plot_angular_momentum(path: Path | None = None) -> None:
    """
    Plot angular momentum evolution over time for an orbit in a logarithmic potential.

    Generates a two-panel plot showing the total angular momentum L(t) over time.
    The total angular momentum includes both the z-component (Lz) and contributions
    from meridional motion. In axisymmetric potentials, only Lz is conserved,
    while L varies due to coupling between radial and vertical motion.

    Parameters
    ----------
    path : Path | None, optional
        File path to save the figure. If None, displays the plot interactively.

    Notes
    -----
    The total angular momentum is calculated as:
    L = √[(z*pR - R*pz)^2 + (R^2 + z^2)*Lz^2/R^2]

    Fixed parameters used:
    - Logarithmic potential with q=0.9, v0=1.0
    - Specific angular momentum Lz=0.2
    - Total energy E=-0.8
    - Initial conditions: R0=0.35, pR0=0.1
    - Integration time: t_max=300.0, dt=0.01

    Panel layout:
    - Left panel: Full time evolution (t=0 to 20)
    - Right panel: Detailed view of late-time behavior (t=200 to 220)

    The variation in L reflects the non-integrability of motion in the
    meridional plane, even though Lz remains constant.
    """
    fig: Figure
    axs: AxesArray
    fig, axs = plt.subplots(
        1,
        2,
        figsize=(6, 5),
        sharex=False,
    )

    fig.subplots_adjust(
        left=0.1,
        right=0.98,
        bottom=0.1,
        top=0.92,
        wspace=0.15,
        hspace=0.05,
    )

    model = LogarithmicPotential(LogarithmicPotentialParams(q=0.9, v0=1.0, Rc=0.0), Lz=0.2)
    t, R, z, pR, pz = integrate_orbit(  # pylint: disable=unused-variable
        model,
        E=-0.8,
        R0=0.35,
        pR0=0.1,
        t_max=300.0,
        dt=0.01,
    )
    L = np.sqrt((z * pR - R * pz) ** 2 + (R**2 + z**2) * model.Lz**2 / R**2)

    # left panel
    axs[0].tick_params(which="minor", length=3, color="gray", direction="in")
    axs[0].tick_params(which="major", length=6, direction="in")
    axs[0].tick_params(top=True, right=True, which="both")

    axs[0].set_xlim(0, 20)
    axs[0].set_ylim(0, 0.3)
    axs[0].set_xlabel(r"$t$")
    axs[0].set_ylabel(r"$L$")
    axs[0].plot(t, L, "k-", label="$L$")

    # right panel
    axs[1].tick_params(which="minor", length=3, color="gray", direction="in")
    axs[1].tick_params(which="major", length=6, direction="in")
    axs[1].tick_params(top=True, right=True, which="both")

    axs[1].set_xlim(200, 220)
    axs[1].set_ylim(0, 0.3)
    axs[1].set_xlabel(r"$t$")
    axs[1].plot(t, L, "k-", label="$L$")

    if path:
        fig.savefig(
            path,
            dpi=150,
            bbox_inches="tight",
            facecolor="white",
            edgecolor="none",
        )
    else:
        plt.show()
