"""Module for calculating properties of an exponential disk."""

import logging
import warnings

import numpy as np
from scipy.integrate import quad
from scipy.special import ive, kve  # pylint: disable=no-name-in-module

warnings.filterwarnings("error")
logger = logging.getLogger(__name__)


def disk_potential(R: float) -> float:
    """
    Calculate the gravitational potential of an exponential disk at radius R.

    Uses a robust approach that avoids discontinuities and overflow issues.
    1. Uses exponentially scaled modified Bessel functions to maintain numerical stability.
    2. Falls back to asymptotic approximations for very large R to prevent overflow

    Parameters
    ----------
    R : float
        Radial distance from the center of the disk.

    Returns
    -------
    float
        Gravitational potential at radius R.
    """
    x = R / 2

    try:

        I0_scaled = ive(0, x)
        K0_scaled = kve(0, x)
        I1_scaled = ive(1, x)
        K1_scaled = kve(1, x)

        combination = I0_scaled * K1_scaled - I1_scaled * K0_scaled
        result = -np.pi * R * combination

        return float(result)

    except (Warning, RuntimeError, OverflowError) as e:
        logger.error("Error in disk_potential: %s for R=%f", e, R)
        return float(-np.pi * np.sqrt(np.pi / (4 * x)) * np.exp(-2 * x))


def disk_density(R: float) -> float:
    """Calculate the mass density of an exponential disk at radius R.

    Parameters
    ----------
    R : float
        Radial distance from the center of the disk.

    Returns
    -------
    float
        Mass density at radius R.
    """
    return float(np.exp(-R))


def potential_energy() -> float:
    """
    Calculate the total potential energy of the exponential disk.

    Returns
    -------
    float
        Total potential energy of the disk.
    """

    def integrand(R: float) -> float:
        return np.pi * R * disk_density(R) * disk_potential(R)

    U, _ = quad(integrand, 0, np.inf)
    return float(U)


def circular_velocity(R: float) -> float:
    """
    Calculate the circular velocity at radius R in the exponential disk.

    Parameters
    ----------
    R : float
        Radial distance from the center of the disk.

    Returns
    -------
    float
        Circular velocity at radius R.
    """
    x = R / 2

    try:

        I0_scaled = ive(0, x)
        K0_scaled = kve(0, x)
        I1_scaled = ive(1, x)
        K1_scaled = kve(1, x)

        combination = I0_scaled * K0_scaled - I1_scaled * K1_scaled
        result = 2 * np.sqrt(np.pi) * x * np.sqrt(combination)

        return float(result)

    except (Warning, RuntimeError, OverflowError) as e:
        logger.error("Error in circular_velocity: %s for R=%f", e, R)
        return float(np.sqrt(2 * np.pi * x) * np.exp(-x))


def angular_momentum() -> float:
    """
    Calculate the total angular momentum of the exponential disk.

    Assuming circular orbits.

    Returns
    -------
    float
        Total angular momentum of the disk.
    """

    def integrand(R: float) -> float:
        return 2 * np.pi * R**2 * disk_density(R) * circular_velocity(R)

    L, _ = quad(integrand, 0, np.inf)
    return float(L)


def kinetic_energy() -> float:
    """
    Calculate the total kinetic energy of the exponential disk.

    Assuming circular orbits.

    Returns
    -------
    float
        Total kinetic energy of the disk.
    """

    def integrand(R: float) -> float:
        v_c = circular_velocity(R)
        return np.pi * R * disk_density(R) * v_c**2

    K, _ = quad(integrand, 0, np.inf)
    return float(K)
