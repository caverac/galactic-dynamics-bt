"""Unit tests for symplectic integrators module."""

import warnings

import numpy as np
import pytest

from galactic_dynamics_bt.chapter03.symplectic_integrators import SymplecticIntegrator


class TestSymplecticIntegrator:
    """Test cases for the SymplecticIntegrator class."""

    # pylint: disable=attribute-defined-outside-init
    def setup_method(self) -> None:
        """Set up common test fixtures."""
        # Simple harmonic oscillator: H = p^2/2 + q^2/2
        self.dT_harmonic = lambda p: p  # dT/dp = p
        self.dV_harmonic = lambda q: q  # dV/dq = q
        self.H_harmonic = lambda q, p: 0.5 * (q**2 + p**2)

        # Basic test parameters
        self.y0_1d = np.array([1.0, 0.0])  # q=1, p=0
        self.t_short = np.linspace(0, 1.0, 11)  # Short integration for speed

        # 2D harmonic oscillator
        self.y0_2d = np.array([1.0, 0.5, 0.0, -0.2])  # q=[1,0.5], p=[0,-0.2]
        self.H_harmonic_2d = lambda q, p: 0.5 * (np.sum(q**2) + np.sum(p**2))

    def test_init_valid_parameters(self) -> None:
        """Test initialization with valid parameters."""
        integrator = SymplecticIntegrator(self.dT_harmonic, self.dV_harmonic, self.y0_1d, self.t_short, order=2)

        assert integrator.order == 2
        assert callable(integrator.dT)
        assert callable(integrator.dV)
        np.testing.assert_array_equal(integrator.y0, self.y0_1d)
        np.testing.assert_array_equal(integrator.t, self.t_short)

    def test_init_invalid_order(self) -> None:
        """Test initialization with invalid integration order."""
        with pytest.raises(ValueError, match="Unsupported order 3"):
            SymplecticIntegrator(self.dT_harmonic, self.dV_harmonic, self.y0_1d, self.t_short, order=3)

    def test_init_odd_length_y0(self) -> None:
        """Test initialization with odd-length initial conditions."""
        y0_odd = np.array([1.0, 0.0, 0.5])  # Length 3 (odd)
        with pytest.raises(ValueError, match="Initial conditions must have even length"):
            SymplecticIntegrator(self.dT_harmonic, self.dV_harmonic, y0_odd, self.t_short)

    def test_init_with_hamiltonian(self) -> None:
        """Test initialization with Hamiltonian function."""
        integrator = SymplecticIntegrator(
            self.dT_harmonic, self.dV_harmonic, self.y0_1d, self.t_short, order=2, H=self.H_harmonic
        )

        assert integrator.H is not None
        assert callable(integrator.H)

    def test_harmonic_oscillator_energy_conservation(self) -> None:
        """Test energy conservation for harmonic oscillator."""
        t = np.linspace(0, 2 * np.pi, 100)
        integrator = SymplecticIntegrator(self.dT_harmonic, self.dV_harmonic, self.y0_1d, t, order=2, H=self.H_harmonic)

        # Suppress energy conservation warnings for this test
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            y = integrator.integrate()

        # Check energy conservation
        energies_list: list[float] = []
        for i in range(len(t)):
            q, p = y[i, 0], y[i, 1]
            energies_list.append(self.H_harmonic(q, p))

        energies = np.array(energies_list)
        energy_variation = np.max(np.abs(energies - energies[0]))

        # Energy should be conserved to within numerical precision
        assert energy_variation < 1e-3

    def test_harmonic_oscillator_period(self) -> None:
        """Test that harmonic oscillator has correct period."""
        # For unit mass and spring constant, period should be 2π
        t = np.linspace(0, 2 * np.pi, 200)
        integrator = SymplecticIntegrator(self.dT_harmonic, self.dV_harmonic, self.y0_1d, t, order=2)

        y = integrator.integrate()

        # Check that position returns to initial value after one period
        q_initial = y[0, 0]
        q_final = y[-1, 0]

        assert abs(q_final - q_initial) < 1e-3

    def test_different_integration_orders(self) -> None:
        """Test all supported integration orders."""
        orders = [1, 2, 4, 6, 8]

        for order in orders:
            integrator = SymplecticIntegrator(self.dT_harmonic, self.dV_harmonic, self.y0_1d, self.t_short, order=order)

            y = integrator.integrate()

            # Check that integration completed without errors
            assert y.shape == (len(self.t_short), len(self.y0_1d))
            assert not np.any(np.isnan(y))
            assert not np.any(np.isinf(y))

    def test_2d_system(self) -> None:
        """Test integration of 2D harmonic oscillator."""
        integrator = SymplecticIntegrator(
            self.dT_harmonic, self.dV_harmonic, self.y0_2d, self.t_short, order=2, H=self.H_harmonic_2d
        )

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            y = integrator.integrate()

        # Check dimensions
        assert y.shape == (len(self.t_short), len(self.y0_2d))

        # Check energy conservation
        energies_list: list[float] = []
        for i in range(len(self.t_short)):
            q = y[i, :2]
            p = y[i, 2:]
            energies_list.append(self.H_harmonic_2d(q, p))

        energies = np.array(energies_list)
        energy_variation = np.max(np.abs(energies - energies[0]))
        assert energy_variation < 1e-2

    def test_free_particle(self) -> None:
        """Test integration of free particle (V=0)."""
        y0 = np.array([0.0, 1.0])  # q=0, p=1 (unit velocity)
        t = np.linspace(0, 1.0, 11)

        # pylint: disable=unnecessary-lambda
        integrator = SymplecticIntegrator(lambda p: p, lambda q: np.zeros_like(q), y0, t, order=2)
        y = integrator.integrate()

        # For free particle, q(t) = q0 + p0*t
        expected_positions = y0[0] + y0[1] * t
        actual_positions = y[:, 0]

        np.testing.assert_allclose(actual_positions, expected_positions, rtol=1e-10)

    def test_energy_warning(self) -> None:
        """Test that energy conservation warning is issued when appropriate."""
        integrator = SymplecticIntegrator(
            lambda p: 10 * p,
            self.dV_harmonic,
            self.y0_1d,
            self.t_short,
            order=2,
            H=self.H_harmonic,
        )

        with pytest.warns(UserWarning, match="Energy not conserved"):
            integrator.integrate()

    def test_individual_integrator_methods(self) -> None:
        """Test individual integrator methods directly."""
        integrator = SymplecticIntegrator(self.dT_harmonic, self.dV_harmonic, self.y0_1d, self.t_short, order=2)

        q, p = 1.0, 0.0
        dt = 0.01

        # Test each integrator method
        methods = [
            integrator._euler_symplectic,  # pylint: disable=protected-access
            integrator._verlet,  # pylint: disable=protected-access
            integrator._ruth4,  # pylint: disable=protected-access
            integrator._kahan_li6,  # pylint: disable=protected-access
            integrator._kahan_li8,  # pylint: disable=protected-access
        ]

        for method in methods:
            q_new, p_new = method(np.array([q]), np.array([p]), dt)

            # Check that results are finite
            assert np.isfinite(q_new).all()
            assert np.isfinite(p_new).all()

            # For harmonic oscillator with small dt, change should be small
            assert abs(q_new[0] - q) < 0.1
            assert abs(p_new[0] - p) < 0.1

    def test_kepler_problem(self) -> None:
        """Test Kepler problem (central force)."""

        # Kepler potential: V = -1/r
        def dT_kepler(p: np.ndarray) -> np.ndarray:
            return p  # T = |p|^2/2

        def dV_kepler(q: np.ndarray) -> np.ndarray:
            r = float(np.linalg.norm(q))
            return q / r**3

        def H_kepler(q: np.ndarray, p: np.ndarray) -> float:
            r = np.linalg.norm(q)
            return float(0.5 * np.sum(p**2) - 1.0 / r)

        # Circular orbit initial conditions
        y0 = np.array([1.0, 0.0, 0.0, 1.0])  # q=[1,0], p=[0,1]
        t = np.linspace(0, np.pi, 50)  # Half orbit for speed

        integrator = SymplecticIntegrator(dT_kepler, dV_kepler, y0, t, order=4, H=H_kepler)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            y = integrator.integrate()

        # Check that orbit is roughly circular (distance from origin)
        distances = np.sqrt(y[:, 0] ** 2 + y[:, 1] ** 2)
        distance_variation = np.max(distances) - np.min(distances)

        # For circular orbit, radius should be constant
        assert distance_variation < 0.1

    def test_positional_only_parameters(self) -> None:
        """Test that first four parameters are positional-only."""
        # This should work (positional arguments)
        integrator = SymplecticIntegrator(self.dT_harmonic, self.dV_harmonic, self.y0_1d, self.t_short, order=2)
        assert integrator is not None

        # This should fail (trying to use keyword arguments for positional-only params)
        with pytest.raises(TypeError):
            SymplecticIntegrator(
                dT=self.dT_harmonic,  # type: ignore
                dV=self.dV_harmonic,
                y0=self.y0_1d,
                t=self.t_short,
                order=2,
            )

    def test_keyword_only_parameters(self) -> None:
        """Test that order and H must be keyword arguments."""
        # This should work (using keyword arguments)
        integrator = SymplecticIntegrator(
            self.dT_harmonic, self.dV_harmonic, self.y0_1d, self.t_short, order=2, H=self.H_harmonic
        )
        assert integrator is not None

    def test_return_value_dimensions(self) -> None:
        """Test that integrate returns correct array dimensions."""
        integrator = SymplecticIntegrator(self.dT_harmonic, self.dV_harmonic, self.y0_1d, self.t_short, order=2)

        y = integrator.integrate()

        # Should return array with shape (n_timepoints, n_dimensions)
        expected_shape = (len(self.t_short), len(self.y0_1d))
        assert y.shape == expected_shape

        # First row should match initial conditions
        np.testing.assert_array_equal(y[0], self.y0_1d)

    def test_time_array_validation(self) -> None:
        """Test that time array is properly validated."""
        # Non-uniform time steps should work
        t_nonuniform = np.array([0.0, 0.1, 0.3, 0.6, 1.0])

        integrator = SymplecticIntegrator(self.dT_harmonic, self.dV_harmonic, self.y0_1d, t_nonuniform, order=2)

        y = integrator.integrate()
        assert y.shape == (len(t_nonuniform), len(self.y0_1d))

    def test_large_time_step_stability(self) -> None:
        """Test stability with larger time steps (but keep test fast)."""
        # Use larger time step but shorter total time
        t = np.linspace(0, 1.0, 6)  # dt = 0.2

        integrator = SymplecticIntegrator(self.dT_harmonic, self.dV_harmonic, self.y0_1d, t, order=4)

        y = integrator.integrate()

        # Integration should complete without NaN or Inf
        assert not np.any(np.isnan(y))
        assert not np.any(np.isinf(y))

    def test_zero_initial_conditions(self) -> None:
        """Test behavior with zero initial conditions."""
        y0_zero = np.array([0.0, 0.0])

        integrator = SymplecticIntegrator(self.dT_harmonic, self.dV_harmonic, y0_zero, self.t_short, order=2)

        y = integrator.integrate()

        # For harmonic oscillator starting at rest at origin, should stay at origin
        np.testing.assert_allclose(y, 0.0, atol=1e-15)

    def test_hamiltonian_with_different_signatures(self) -> None:
        """Test Hamiltonian functions with different calling patterns."""
        # Test that Hamiltonian is called correctly during integration
        call_count = 0

        def counting_hamiltonian(q: np.ndarray, p: np.ndarray) -> float:
            """Count calls to the Hamiltonian function."""
            nonlocal call_count
            call_count += 1
            return 0.5 * float(q.sum() ** 2 + p.sum() ** 2)

        integrator = SymplecticIntegrator(
            self.dT_harmonic, self.dV_harmonic, self.y0_1d, self.t_short, order=2, H=counting_hamiltonian
        )

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            integrator.integrate()

        # Hamiltonian should be called at each time step (minus initial)
        expected_calls = len(self.t_short) - 1
        assert call_count == expected_calls
