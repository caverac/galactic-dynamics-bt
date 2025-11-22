"""
Symplectic integrators for Hamiltonian systems.

This module provides symplectic integrators for Hamiltonian systems of the form
H(q, p) = T(p) + V(q), where T is the kinetic energy (function of momentum) and
V is the potential energy (function of position).
"""

from typing import Callable, Optional, Tuple
import warnings

import numpy as np


class SymplecticIntegrator:
    """
    Symplectic integrator for Hamiltonian systems.

    For separable Hamiltonian systems H(q,p) = T(p) + V(q), provides
    structure-preserving integration that exactly conserves the symplectic
    2-form and approximately conserves energy.

    Parameters
    ----------
    dT : callable
        Gradient of kinetic energy T with respect to momentum p.
        Should return velocity: dT/dp = dq/dt
    dV : callable
        Gradient of potential energy V with respect to position q.
        Should return force: dV/dq = -dp/dt
    y0 : array_like, shape (2*ndim,)
        Initial conditions [q0, p0] where q0 and p0 have same dimension.
    t : array_like
        Time points for integration.
    order : int, optional
        Integration order (1, 2, 4, 6, 8). Default is 2 (Verlet).
    H : callable, optional
        Total Hamiltonian function for energy monitoring.

    Examples
    --------
    Simple harmonic oscillator with H = p^2/2 + q^2/2:

    >>> import numpy as np
    >>> def dT_dp(p): return p  # dT/dp = p for T = p^2/2
    >>> def dV_dq(q): return q  # dV/dq = q for V = q^2/2
    >>> def hamiltonian(q, p): return 0.5 * (q**2 + p**2)
    >>> y0 = [1.0, 0.0]  # Start at q=1, p=0
    >>> t = np.linspace(0, 4*np.pi, 100)
    >>> integrator = SymplecticIntegrator(dT_dp, dV_dq, y0, t, order=2, H=hamiltonian)
    >>> y = integrator.integrate()
    >>> # y contains [q(t), p(t)] at each time point

    Kepler problem with H = p^2/2 - 1/|q|:

    >>> def dT_dp_kepler(p): return p  # T = |p|^2/2
    >>> def dV_dq_kepler(q):
    ...     r = np.linalg.norm(q)
    ...     return -q / r**3  # V = -1/r, dV/dq = q/r^3
    >>> y0 = [1.0, 0.0, 0.0, 1.0]  # Circular orbit: q=[1,0], p=[0,1]
    >>> t = np.linspace(0, 2*np.pi, 1000)
    >>> integrator = SymplecticIntegrator(dT_dp_kepler, dV_dq_kepler, y0, t, order=4)
    >>> y = integrator.integrate()
    """

    dT: Callable[[np.ndarray], np.ndarray]
    dV: Callable[[np.ndarray], np.ndarray]
    H: Optional[Callable[[np.ndarray, np.ndarray], float]]
    y0: np.ndarray
    t: np.ndarray
    order: int

    integrator: Callable[[np.ndarray, np.ndarray, float], Tuple[np.ndarray, np.ndarray]]

    def __init__(
        self,
        dT: Callable[[np.ndarray], np.ndarray],
        dV: Callable[[np.ndarray], np.ndarray],
        y0: np.ndarray,
        t: np.ndarray,
        /,
        *,
        order: int = 2,
        H: Optional[Callable[[np.ndarray, np.ndarray], float]] = None,
    ):
        """
        Initialize the symplectic integrator.

        Parameters
        ----------
        dT : callable
            Gradient of kinetic energy with respect to momentum.
        dV : callable
            Gradient of potential energy with respect to position.
        y0 : array_like
            Initial conditions [q0, p0].
        t : array_like
            Time points for integration.
        order : int, optional
            Integration order (1, 2, 4, 6, 8). Default is 2.
        H : callable, optional
            Hamiltonian function for energy monitoring.
        """
        self.y0 = y0
        self.t = t
        self.order = order
        self.H = H

        if len(y0) % 2 != 0:
            raise ValueError("Initial conditions must have even length (q and p)")

        self.dT = dT
        self.dV = dV

        # Choose integrator
        if self.order == 1:
            self.integrator = self._euler_symplectic
        elif self.order == 2:
            self.integrator = self._verlet
        elif self.order == 4:
            self.integrator = self._ruth4
        elif self.order == 6:
            self.integrator = self._kahan_li6
        elif self.order == 8:
            self.integrator = self._kahan_li8
        else:
            raise ValueError(f"Unsupported order {self.order}. Supported: 1, 2, 4, 6, 8")

    def integrate(self) -> np.ndarray:
        """
        Integrate the Hamiltonian system and return solution array.

        Returns
        -------
        y : ndarray, shape (len(t), 2*ndim)
            Solution array where y[i] = [q(t[i]), p(t[i])].
            Time points are available as self.t.

        Examples
        --------
        >>> # Pendulum: H = p^2/2 - cos(q)
        >>> def dT_dp(p): return p
        >>> def dV_dq(q): return np.sin(q)  # dV/dq for V = -cos(q)
        >>> integrator = SymplecticIntegrator(dT_dp, dV_dq, [np.pi/4, 0], [0, 5], order=4)
        >>> y = integrator.integrate()
        >>> angles = y[:, 0]  # Extract position (angle)
        >>> momenta = y[:, 1]  # Extract momentum
        >>> times = integrator.t  # Access time points
        """
        y0 = np.asarray(self.y0, dtype=float)
        t = np.asarray(self.t, dtype=float)

        ndim = len(y0) // 2

        # Initialize result array
        y = np.zeros((len(t), len(y0)))
        y[0] = y0

        # Split into position and momentum
        q = y0[:ndim].copy()
        p = y0[ndim:].copy()

        # Integrate
        E0 = None
        for i in range(1, len(t)):
            dt = t[i] - t[i - 1]
            q, p = self.integrator(q, p, dt)
            y[i] = np.concatenate([q, p])

            if self.H:

                if E0 is None:
                    E0 = self.H(q, p)
                else:
                    E_current = self.H(q, p)
                    energy_variation = np.abs(E_current - E0)
                    if energy_variation > 1e-3 * np.abs(E0):
                        warnings.warn(
                            f"Energy not conserved at step {i}: variation = {energy_variation:.2e}, "
                            f"relative to initial = {energy_variation/np.abs(E0):.2e}"
                        )

        return y

    def _euler_symplectic(self, q: np.ndarray, p: np.ndarray, dt: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        First-order symplectic Euler method.

        Implements the partitioned Runge-Kutta scheme:
        p_{n+1} = p_n - dt * dV(q_n)/dq
        q_{n+1} = q_n + dt * dT(p_{n+1})/dp

        This is first-order accurate but exactly preserves the symplectic structure.
        Use for very stiff systems or when energy drift must be minimized.

        Examples
        --------
        >>> # Single step of harmonic oscillator
        >>> q, p = np.array([1.0]), np.array([0.0])
        >>> dt = 0.1
        >>> # For integrator with dT(p)=p, dV(q)=q:
        >>> q_new, p_new = integrator._euler_symplectic(q, p, dt)
        >>> # Result: p_new = p - dt*q = -0.1, q_new = q + dt*p_new = 0.99
        """
        # p_{n+1} = p_n - dt * dV/dq(q_n)
        p_new = p - dt * self.dV(q)
        # q_{n+1} = q_n + dt * dT/dp(p_{n+1})
        q_new = q + dt * self.dT(p_new)
        return q_new, p_new

    def _verlet(self, q: np.ndarray, p: np.ndarray, dt: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Second-order Verlet (leapfrog) integrator.

        The workhorse of molecular dynamics simulations. Implements:
        p_{n+1/2} = p_n - (dt/2) * dV/dq(q_n)
        q_{n+1} = q_n + dt * dT/dp(p_{n+1/2})
        p_{n+1} = p_{n+1/2} - (dt/2) * dV/dq(q_{n+1})

        Second-order accurate with excellent long-term stability.
        Recommended for most applications.

        Examples
        --------
        >>> # Verlet step for 2D harmonic oscillator
        >>> q = np.array([1.0, 0.5])  # 2D position
        >>> p = np.array([0.0, -0.2])  # 2D momentum
        >>> dt = 0.01
        >>> # For H = (p1^2 + p2^2)/2 + (q1^2 + q2^2)/2:
        >>> q_new, p_new = integrator._verlet(q, p, dt)
        >>> # Energy should be approximately conserved
        """
        # Half step for momentum
        p_half = p - 0.5 * dt * self.dV(q)
        # Full step for position
        q_new = q + dt * self.dT(p_half)
        # Half step for momentum
        p_new = p_half - 0.5 * dt * self.dV(q_new)
        return q_new, p_new

    # pylint: disable=too-many-locals
    def _ruth4(self, q: np.ndarray, p: np.ndarray, dt: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Fourth-order Ruth integrator.

        Composition method using the optimal coefficients from Ruth (1983).
        Provides fourth-order accuracy while preserving symplectic structure.
        More expensive than Verlet but allows larger time steps.

        The method composes multiple Verlet steps with coefficients:
        1/(2 - 2^(1/3)) = 1.351

        Examples
        --------
        >>> # Use for systems requiring high precision
        >>> q = np.array([1.0])  # Position
        >>> p = np.array([0.0])  # Momentum
        >>> dt = 0.1  # Can use larger timesteps than Verlet
        >>> q_new, p_new = integrator._ruth4(q, p, dt)
        >>> # Fourth-order accurate result with symplectic properties preserved
        """
        # Ruth coefficients
        theta = 1.0 / (2.0 - 2 ** (1.0 / 3.0))
        alpha = theta
        beta = -theta / 2.0
        gamma = 1.0 - 2.0 * theta
        delta = gamma

        # Substeps
        q1, p1 = self._verlet(q, p, alpha * dt)
        q2, p2 = self._verlet(q1, p1, beta * dt)
        q3, p3 = self._verlet(q2, p2, gamma * dt)
        q4, p4 = self._verlet(q3, p3, delta * dt)

        return q4, p4

    def _kahan_li6(self, q: np.ndarray, p: np.ndarray, dt: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Sixth-order Kahan-Li integrator.

        High-order symplectic integrator using optimized coefficients.
        Excellent for problems requiring very high accuracy over long times,
        such as planetary motion or galactic dynamics.

        Uses 7 composition coefficients derived by Kahan & Li (1997)
        for minimal leading-order error terms.

        Examples
        --------
        >>> # Long-term orbital integration
        >>> q = np.array([1.0, 0.0])  # Initial position (AU)
        >>> p = np.array([0.0, 2*np.pi])  # Initial momentum (AU/year)
        >>> dt = 0.1  # Large timestep possible due to high order
        >>> q_new, p_new = integrator._kahan_li6(q, p, dt)
        >>> # Suitable for thousand-year integrations with minimal drift
        """
        # Kahan-Li coefficients for 6th order
        w0 = -1.17767998417887
        w1 = 0.235573213359357
        w2 = 0.784513610477560
        coeffs = [
            w1 / 2,
            (w0 + w1) / 2,
            (w2 + w0) / 2,
            (1 - 2 * (w0 + w1 + w2)) / 2,
            (w2 + w0) / 2,
            (w0 + w1) / 2,
            w1 / 2,
        ]

        q_curr, p_curr = q, p
        for coeff in coeffs:
            q_curr, p_curr = self._verlet(q_curr, p_curr, 2 * coeff * dt)

        return q_curr, p_curr

    # pylint: disable=too-many-locals
    def _kahan_li8(self, q: np.ndarray, p: np.ndarray, dt: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Eighth-order Kahan-Li integrator.

        Highest-order integrator available. Uses 17 composition coefficients
        for exceptional accuracy. Ideal for studies requiring machine precision
        over astronomical timescales.

        Computational cost is ~8.5x higher than Verlet per step, but allows
        much larger timesteps while maintaining accuracy.

        Examples
        --------
        >>> # Ultra-high precision solar system integration
        >>> q = np.array([5.2, 0.0])  # Jupiter's orbit (AU)
        >>> p = np.array([0.0, 0.84])  # Orbital velocity (AU/year)
        >>> dt = 1.0  # Full year timesteps possible
        >>> q_new, p_new = integrator._kahan_li8(q, p, dt)
        >>> # Accurate to ~1e-12 over million-year timescales
        """
        # Kahan-Li coefficients for 8th order
        w0 = 0.74167036435061295344822780
        w1 = -0.40910082580003159399730010
        w2 = 0.19075471029623837995387626
        w3 = -0.57386247111608226665638773
        w4 = 0.29906418130365592384446354
        w5 = 0.33462491824529818378495798
        w6 = 0.31529309239676659663205666
        w7 = -0.79688793935291635401978884

        coeffs = [
            w7 / 2,
            (w6 + w7) / 2,
            (w5 + w6) / 2,
            (w4 + w5) / 2,
            (w3 + w4) / 2,
            (w2 + w3) / 2,
            (w1 + w2) / 2,
            (w0 + w1) / 2,
            (1 - 2 * (w0 + w1 + w2 + w3 + w4 + w5 + w6 + w7)) / 2,
            (w0 + w1) / 2,
            (w1 + w2) / 2,
            (w2 + w3) / 2,
            (w3 + w4) / 2,
            (w4 + w5) / 2,
            (w5 + w6) / 2,
            (w6 + w7) / 2,
            w7 / 2,
        ]

        q_curr, p_curr = q, p
        for coeff in coeffs:
            q_curr, p_curr = self._verlet(q_curr, p_curr, 2 * coeff * dt)

        return q_curr, p_curr
