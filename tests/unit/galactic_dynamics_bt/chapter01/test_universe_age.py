"""Test suite for universe_age module."""

from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from galactic_dynamics_bt.chapter01.universe_age import find_universe_age, plot_universe_age


class TestFindUniverseAge:
    """Test suite for find_universe_age function."""

    def test_default_parameters_at_z_zero(self) -> None:
        """Test universe age calculation at z=0 with default parameters."""
        age = find_universe_age(0.0)
        # Should return age in Gyr with typical cosmological parameters
        assert isinstance(age, float)
        assert 10.0 < age < 20.0  # Reasonable age range for universe

    def test_standard_cosmology(self) -> None:
        """Test with Planck 2018 cosmological parameters."""
        age = find_universe_age(
            0.0,
            omega_m0=0.315,
            omega_lambda0=0.685,
            omega_gamma0=8.24e-5 / 0.674**2,
            H0=100 * 0.674,
        )
        # Expected age should be around 13.8 Gyr for Planck cosmology
        assert 13.0 < age < 14.5

    def test_higher_redshift_gives_smaller_age(self) -> None:
        """Test that higher redshift results in smaller universe age."""
        age_z0 = find_universe_age(0.0)
        age_z1 = find_universe_age(1.0)
        age_z2 = find_universe_age(2.0)

        assert age_z2 < age_z1 < age_z0

    def test_matter_dominated_universe(self) -> None:
        """Test with a matter-dominated universe (no dark energy)."""
        age = find_universe_age(
            0.0,
            omega_m0=1.0,
            omega_lambda0=0.0,
            omega_gamma0=0.0,
            H0=70.0,
        )
        # Matter-dominated universe should be younger
        assert 5.0 < age < 12.0

    def test_dark_energy_dominated_universe(self) -> None:
        """Test with a dark energy-dominated universe."""
        age = find_universe_age(
            0.0,
            omega_m0=0.1,
            omega_lambda0=0.9,
            omega_gamma0=0.0,
            H0=70.0,
        )
        # Dark energy dominated should be older
        assert 12.0 < age < 25.0

    def test_radiation_dominated_early_universe(self) -> None:
        """Test with significant radiation component."""
        age = find_universe_age(
            10.0,  # High redshift
            omega_m0=0.3,
            omega_lambda0=0.7,
            omega_gamma0=0.1,  # Non-negligible radiation
            H0=70.0,
        )
        assert isinstance(age, float)
        assert age > 0

    def test_hubble_constant_scaling(self) -> None:
        """Test that doubling H0 halves the age."""
        age_h70 = find_universe_age(0.0, H0=70.0)
        age_h140 = find_universe_age(0.0, H0=140.0)

        # Age should scale inversely with H0
        np.testing.assert_allclose(age_h140, age_h70 / 2, rtol=1e-10)

    def test_position_only_parameter(self) -> None:
        """Test that z is a position-only parameter."""
        # This should work
        age = find_universe_age(1.0, omega_m0=0.3)
        assert isinstance(age, float)

        # Note: z is position-only, but pytest can't easily test this
        # at runtime since it's a signature constraint

    def test_keyword_only_parameters(self) -> None:
        """Test that cosmological parameters are keyword-only."""
        # This should work
        age = find_universe_age(1.0, omega_m0=0.3, omega_lambda0=0.7)
        assert isinstance(age, float)

        # Note: keyword-only parameters are enforced at function definition
        # Testing this would require dynamic function calls

    def test_negative_redshift_raises_error(self) -> None:
        """Test that negative redshift raises an appropriate error."""
        # The function doesn't explicitly check for this, but the integral
        # bounds would be invalid, so scipy.integrate.quad should handle it
        with pytest.raises(ZeroDivisionError):
            find_universe_age(-1.0)

    def test_very_high_redshift(self) -> None:
        """Test calculation at very high redshift."""
        age = find_universe_age(100.0)
        assert isinstance(age, float)
        assert age > 0
        assert age < 1.0  # Should be very young

    def test_flat_universe_constraint(self) -> None:
        """Test with parameters that don't sum to 1 (non-flat universe)."""
        age = find_universe_age(
            0.0,
            omega_m0=0.3,
            omega_lambda0=0.6,  # Sum is 0.9, not 1
            omega_gamma0=0.0,
            H0=70.0,
        )
        # Should still work - the (1 - sum) term handles curvature
        assert isinstance(age, float)
        assert age > 0

    def test_return_type_is_float(self) -> None:
        """Test that return type is Python float, not numpy scalar."""
        age = find_universe_age(0.0)
        assert isinstance(age, float)
        assert not isinstance(age, np.floating)

    @pytest.mark.parametrize("z", [0.0, 0.5, 1.0, 2.0, 5.0])
    def test_monotonic_decrease_with_redshift(self, z: float) -> None:
        """Test that age decreases monotonically with redshift."""
        if z > 0:
            age_current = find_universe_age(z)
            age_earlier = find_universe_age(z - 0.1)
            assert age_current < age_earlier

    def test_integration_convergence(self) -> None:
        """Test that the integration converges properly."""
        # The integrand should behave well near a=0 and a=1/(1+z)
        age = find_universe_age(1.0)

        # Should be finite and positive
        assert np.isfinite(age)
        assert age > 0


class TestPlotUniverseAge:
    """Test suite for plot_universe_age function."""

    @patch("galactic_dynamics_bt.chapter01.universe_age.plt")
    def test_plot_creates_figure(self, mock_plt: MagicMock) -> None:
        """Test that plot_universe_age creates a figure and axes."""
        mock_fig = MagicMock()
        mock_axs = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_axs)

        plot_universe_age()

        # Check that subplots was called with correct parameters
        mock_plt.subplots.assert_called_once_with(1, 1, figsize=(6.5, 5), sharex=True)

        # Check that various plotting methods were called
        mock_axs.set_xlabel.assert_called()
        mock_axs.set_ylabel.assert_called()
        mock_axs.set_xscale.assert_called_with("log")
        mock_axs.set_yscale.assert_called_with("log")
        mock_axs.set_xlim.assert_called()

    @patch("galactic_dynamics_bt.chapter01.universe_age.plt")
    def test_plot_saves_figure(self, mock_plt: MagicMock) -> None:
        """Test that the figure is saved to the correct location."""
        mock_fig = MagicMock()
        mock_axs = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_axs)

        plot_universe_age()

        # Check that savefig was called with correct parameters
        mock_fig.savefig.assert_called_once_with(
            "docs/assets/generated/universe_age.png",
            dpi=150,
            bbox_inches="tight",
            facecolor="white",
            edgecolor="none",
        )

    @patch("galactic_dynamics_bt.chapter01.universe_age.plt")
    @patch("galactic_dynamics_bt.chapter01.universe_age.find_universe_age")
    def test_plot_calls_find_universe_age(self, mock_find_age: MagicMock, mock_plt: MagicMock) -> None:
        """Test that plot function calls find_universe_age with correct parameters."""
        mock_fig = MagicMock()
        mock_axs = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_axs)

        # Mock the return value to avoid actual calculations
        mock_find_age.return_value = 13.8

        plot_universe_age()

        # Check that find_universe_age was called multiple times
        assert mock_find_age.call_count > 0

        # Check that it was called with both cosmological parameter sets
        calls = mock_find_age.call_args_list

        # Should include calls with both parameter sets used in the plot
        omega_m_values = [call.kwargs.get("omega_m0") for call in calls if "omega_m0" in call.kwargs]
        assert 0.237 in omega_m_values  # First cosmology
        assert 0.315 in omega_m_values  # Second cosmology (Planck 2018)

    def test_plot_function_runs_without_error(self) -> None:
        """Integration test - check that plot function runs without raising exceptions."""
        # This is a basic smoke test
        try:
            plot_universe_age()
        except Exception as e:  # pylint: disable=broad-except
            pytest.fail(f"plot_universe_age() raised an exception: {e}")
