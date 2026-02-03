"""Unit tests for the exponential_disk module."""

# pylint: disable=import-outside-toplevel
from typing import Callable
import unittest
from unittest.mock import Mock, patch

import numpy as np

from galactic_dynamics_bt.chapter02.exponential_disk import (
    angular_momentum,
    circular_velocity,
    disk_density,
    disk_potential,
    kinetic_energy,
    potential_energy,
)


class TestDiskDensity(unittest.TestCase):
    """Test the disk_density function."""

    def test_disk_density_at_center(self) -> None:
        """Test disk density at the center (R=0)."""
        result = disk_density(0.0)
        expected = 1.0  # exp(-0) = 1
        self.assertAlmostEqual(result, expected, places=10)

    def test_disk_density_exponential_decay(self) -> None:
        """Test that disk density follows exponential decay."""
        # Test at R = 1
        result = disk_density(1.0)
        expected = np.exp(-1.0)
        self.assertAlmostEqual(result, expected, places=10)

        # Test at R = 2
        result = disk_density(2.0)
        expected = np.exp(-2.0)
        self.assertAlmostEqual(result, expected, places=10)

    def test_disk_density_large_radius(self) -> None:
        """Test disk density at large radius."""
        result = disk_density(10.0)
        expected = np.exp(-10.0)
        self.assertAlmostEqual(result, expected, places=10)
        self.assertLess(result, 1e-4)  # Should be very small

    def test_disk_density_return_type(self) -> None:
        """Test that disk_density returns a float."""
        result = disk_density(1.0)
        self.assertIsInstance(result, float)


class TestDiskPotential(unittest.TestCase):
    """Test the disk_potential function."""

    def test_disk_potential_small_radius(self) -> None:
        """Test disk potential at small radius."""
        result = disk_potential(0.1)
        self.assertIsInstance(result, float)
        self.assertTrue(np.isfinite(result))
        self.assertLess(result, 0)  # Potential should be negative

    def test_disk_potential_moderate_radius(self) -> None:
        """Test disk potential at moderate radius."""
        result = disk_potential(5.0)
        self.assertIsInstance(result, float)
        self.assertTrue(np.isfinite(result))
        self.assertLess(result, 0)  # Potential should be negative

    def test_disk_potential_large_radius(self) -> None:
        """Test disk potential at large radius."""
        result = disk_potential(100.0)
        self.assertIsInstance(result, float)
        self.assertTrue(np.isfinite(result))
        self.assertLess(result, 0)  # Potential should be negative

    def test_disk_potential_very_large_radius(self) -> None:
        """Test disk potential at very large radius (asymptotic regime)."""
        result = disk_potential(1000.0)
        self.assertIsInstance(result, float)
        self.assertTrue(np.isfinite(result))
        self.assertLess(result, 0)  # Potential should be negative

    def test_disk_potential_continuity(self) -> None:
        """Test that disk potential is continuous across different regimes."""
        # Test continuity around the transition region
        R_values = np.linspace(90, 110, 21)
        potentials = [disk_potential(R) for R in R_values]

        # Check that all values are finite
        for pot in potentials:
            self.assertTrue(np.isfinite(pot))

        # Check for smoothness - no large jumps
        for i in range(1, len(potentials)):
            relative_change = abs((potentials[i] - potentials[i - 1]) / potentials[i - 1])
            self.assertLess(relative_change, 0.1)  # Less than 10% change per step

    def test_disk_potential_monotonic_behavior(self) -> None:
        """Test that disk potential behaves monotonically in asymptotic regime."""
        # At large radii, potential should approach zero from below
        R_large = [200, 400, 800]
        potentials = [disk_potential(R) for R in R_large]

        # Each should be closer to zero than the previous
        for i in range(1, len(potentials)):
            self.assertGreater(potentials[i], potentials[i - 1])  # More negative means smaller
            self.assertLess(potentials[i], 0)  # Still negative

    @patch("galactic_dynamics_bt.chapter02.exponential_disk.logger")
    def test_disk_potential_error_handling(self, mock_logger: Mock) -> None:
        """Test error handling in disk_potential."""
        # This test checks that the function handles numerical issues gracefully
        with patch("galactic_dynamics_bt.chapter02.exponential_disk.ive", side_effect=RuntimeError("Test error")):
            result = disk_potential(50.0)
            self.assertIsInstance(result, float)
            self.assertTrue(np.isfinite(result))
            mock_logger.error.assert_called()


class TestCircularVelocity(unittest.TestCase):
    """Test the circular_velocity function."""

    def test_circular_velocity_small_radius(self) -> None:
        """Test circular velocity at small radius."""
        result = circular_velocity(0.1)
        self.assertIsInstance(result, float)
        self.assertTrue(np.isfinite(result))
        self.assertGreater(result, 0)  # Velocity should be positive

    def test_circular_velocity_moderate_radius(self) -> None:
        """Test circular velocity at moderate radius."""
        result = circular_velocity(5.0)
        self.assertIsInstance(result, float)
        self.assertTrue(np.isfinite(result))
        self.assertGreater(result, 0)  # Velocity should be positive

    def test_circular_velocity_large_radius(self) -> None:
        """Test circular velocity at large radius."""
        result = circular_velocity(100.0)
        self.assertIsInstance(result, float)
        self.assertTrue(np.isfinite(result))
        self.assertGreater(result, 0)  # Velocity should be positive

    def test_circular_velocity_asymptotic_behavior(self) -> None:
        """Test that circular velocity approaches zero at large radii."""
        v1 = circular_velocity(100.0)
        v2 = circular_velocity(200.0)
        v3 = circular_velocity(400.0)

        # Velocity should decrease with radius in the asymptotic regime
        self.assertGreater(v1, v2)
        self.assertGreater(v2, v3)

    @patch("galactic_dynamics_bt.chapter02.exponential_disk.logger")
    def test_circular_velocity_error_handling(self, mock_logger: Mock) -> None:
        """Test error handling in circular_velocity."""
        with patch("galactic_dynamics_bt.chapter02.exponential_disk.ive", side_effect=RuntimeError("Test error")):
            result = circular_velocity(50.0)
            self.assertIsInstance(result, float)
            self.assertTrue(np.isfinite(result))
            mock_logger.error.assert_called()


class TestIntegratedQuantities(unittest.TestCase):
    """Test the integrated quantities (potential energy, angular momentum, kinetic energy)."""

    @patch("galactic_dynamics_bt.chapter02.exponential_disk.quad")
    def test_potential_energy_integration(self, mock_quad: Mock) -> None:
        """Test potential energy calculation with mocked integration."""
        # Mock the quad function to return a known value
        mock_quad.return_value = (-10.5, 1e-10)

        result = potential_energy()
        self.assertIsInstance(result, float)
        self.assertEqual(result, -10.5)

        # Verify quad was called with correct arguments
        mock_quad.assert_called_once()
        args, _ = mock_quad.call_args
        self.assertEqual(len(args), 3)  # integrand, 0, np.inf
        self.assertEqual(args[1], 0)
        self.assertEqual(args[2], np.inf)

    @patch("galactic_dynamics_bt.chapter02.exponential_disk.quad")
    def test_angular_momentum_integration(self, mock_quad: Mock) -> None:
        """Test angular momentum calculation with mocked integration."""
        # Mock the quad function to return a known value
        mock_quad.return_value = (15.2, 1e-10)

        result = angular_momentum()
        self.assertIsInstance(result, float)
        self.assertEqual(result, 15.2)

        # Verify quad was called
        mock_quad.assert_called_once()

    @patch("galactic_dynamics_bt.chapter02.exponential_disk.quad")
    def test_kinetic_energy_integration(self, mock_quad: Mock) -> None:
        """Test kinetic energy calculation with mocked integration."""
        # Mock the quad function to return a known value
        mock_quad.return_value = (8.7, 1e-10)

        result = kinetic_energy()
        self.assertIsInstance(result, float)
        self.assertEqual(result, 8.7)

        # Verify quad was called
        mock_quad.assert_called_once()

    def test_potential_energy_real_calculation(self) -> None:
        """Test potential energy with a small integration range for speed."""
        # This is a real calculation but with limited range for testing
        with patch("galactic_dynamics_bt.chapter02.exponential_disk.quad") as mock_quad:
            # Use a small finite upper limit for testing
            def mock_quad_finite(
                integrand: Callable[[float], float], lower: float, upper: float
            ) -> tuple[float, float]:
                """Mock quad function with finite upper limit for testing."""
                from scipy.integrate import quad as real_quad

                if upper == np.inf:
                    # Replace infinity with a reasonable upper limit for testing
                    return real_quad(integrand, lower, 10.0)
                return real_quad(integrand, lower, upper)

            mock_quad.side_effect = mock_quad_finite

            result = potential_energy()
            self.assertIsInstance(result, float)
            self.assertTrue(np.isfinite(result))
            self.assertLess(result, 0)  # Potential energy should be negative

    def test_angular_momentum_real_calculation(self) -> None:
        """Test angular momentum with a small integration range for speed."""
        with patch("galactic_dynamics_bt.chapter02.exponential_disk.quad") as mock_quad:

            def mock_quad_finite(
                integrand: Callable[[float], float], lower: float, upper: float
            ) -> tuple[float, float]:
                from scipy.integrate import quad as real_quad

                if upper == np.inf:
                    return real_quad(integrand, lower, 10.0)
                return real_quad(integrand, lower, upper)

            mock_quad.side_effect = mock_quad_finite

            result = angular_momentum()
            self.assertIsInstance(result, float)
            self.assertTrue(np.isfinite(result))
            self.assertGreater(result, 0)  # Angular momentum should be positive

    def test_kinetic_energy_real_calculation(self) -> None:
        """Test kinetic energy with a small integration range for speed."""
        with patch("galactic_dynamics_bt.chapter02.exponential_disk.quad") as mock_quad:

            def mock_quad_finite(
                integrand: Callable[[float], float], lower: float, upper: float
            ) -> tuple[float, float]:
                from scipy.integrate import quad as real_quad

                if upper == np.inf:
                    return real_quad(integrand, lower, 10.0)
                return real_quad(integrand, lower, upper)

            mock_quad.side_effect = mock_quad_finite

            result = kinetic_energy()
            self.assertIsInstance(result, float)
            self.assertTrue(np.isfinite(result))
            self.assertGreater(result, 0)  # Kinetic energy should be positive


class TestPhysicalConsistency(unittest.TestCase):
    """Test physical consistency of the exponential disk model."""

    def test_virial_theorem_approximation(self) -> None:
        """Test that the virial theorem is approximately satisfied."""
        # For testing purposes, use limited integration range
        with patch("galactic_dynamics_bt.chapter02.exponential_disk.quad") as mock_quad:

            def mock_quad_finite(
                integrand: Callable[[float], float], lower: float, upper: float
            ) -> tuple[float, float]:
                from scipy.integrate import quad as real_quad

                if upper == np.inf:
                    return real_quad(integrand, lower, 10.0)
                return real_quad(integrand, lower, upper)

            mock_quad.side_effect = mock_quad_finite

            K = kinetic_energy()
            U = potential_energy()

            # For a disk in equilibrium, 2K + U should be approximately zero (virial theorem)
            virial_ratio = abs((2 * K + U) / K)
            self.assertLess(virial_ratio, 1.0)  # Should be reasonably close to zero

    def test_mass_conservation(self) -> None:
        """Test that the integrated mass is finite and positive."""
        # Mass = integral of 2*pi * R * rho(R) dR
        with patch("galactic_dynamics_bt.chapter02.exponential_disk.quad") as mock_quad:

            def mock_mass_integrand(R: float) -> float:
                """Integrand for mass calculation."""
                return 2 * np.pi * R * disk_density(R)

            def mock_quad_finite(
                integrand: Callable[[float], float], lower: float, upper: float
            ) -> tuple[float, float]:
                """Mock quad function with finite upper limit for testing."""
                from scipy.integrate import quad as real_quad

                if upper == np.inf:
                    return real_quad(integrand, lower, 10.0)
                return real_quad(integrand, lower, upper)

            mock_quad.side_effect = mock_quad_finite

            from scipy.integrate import quad

            mass, _ = quad(mock_mass_integrand, 0, 10.0)

            self.assertGreater(mass, 0)
            self.assertTrue(np.isfinite(mass))


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and boundary conditions."""

    def test_zero_radius(self) -> None:
        """Test behavior at R = 0."""
        # Density should be 1 at center
        density = disk_density(0.0)
        self.assertEqual(density, 1.0)

        # Potential and velocity should be finite at center
        # Note: potential might be singular at R=0, so we test a very small value
        potential = disk_potential(1e-6)
        velocity = circular_velocity(1e-6)

        self.assertTrue(np.isfinite(potential))
        self.assertTrue(np.isfinite(velocity))

    def test_very_small_radius(self) -> None:
        """Test behavior at very small but non-zero radius."""
        R_small = 1e-10

        density = disk_density(R_small)
        potential = disk_potential(R_small)
        velocity = circular_velocity(R_small)

        self.assertTrue(np.isfinite(density))
        self.assertTrue(np.isfinite(potential))
        self.assertTrue(np.isfinite(velocity))

        self.assertAlmostEqual(density, 1.0, places=8)  # Should be very close to 1

    def test_consistency_across_scales(self) -> None:
        """Test that functions give consistent results across different scales."""
        R_values = [0.01, 0.1, 1.0, 10.0, 100.0]

        for R in R_values:
            density = disk_density(R)
            potential = disk_potential(R)
            velocity = circular_velocity(R)

            # All should be finite
            self.assertTrue(np.isfinite(density), f"Density not finite at R={R}")
            self.assertTrue(np.isfinite(potential), f"Potential not finite at R={R}")
            self.assertTrue(np.isfinite(velocity), f"Velocity not finite at R={R}")

            # Physical constraints
            self.assertGreater(density, 0, f"Density not positive at R={R}")
            self.assertLess(potential, 0, f"Potential not negative at R={R}")
            self.assertGreater(velocity, 0, f"Velocity not positive at R={R}")


if __name__ == "__main__":
    unittest.main()
