"""Unit tests for Sersic profile functionality."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from galactic_dynamics_bt.chapter01.sersic_profile import get_luminosity_sersic, plot_sersic_profile, sersic_profile


class TestSersicProfile:
    """Test cases for Sersic profile functions."""

    def test_sersic_profile_basic_functionality(self) -> None:
        """Test basic functionality of sersic_profile function."""
        # Test with default parameters
        I, dI_dR = sersic_profile(1.0)

        # Should return floats
        assert isinstance(I, float)
        assert isinstance(dI_dR, float)

        # At Re=1, I should equal Ie=1 (definition of effective radius)
        assert np.isclose(I, 1.0, rtol=1e-10)

        # Derivative should be negative (decreasing profile)
        assert dI_dR < 0

    def test_sersic_profile_parameter_variations(self) -> None:
        """Test sersic_profile with different parameters."""
        # Test different Sersic indices
        I_exp, _ = sersic_profile(1.0, m=1)  # Exponential
        I_dev, _ = sersic_profile(1.0, m=4)  # de Vaucouleurs

        # At effective radius, both should equal Ie
        assert np.isclose(I_exp, 1.0, rtol=1e-10)
        assert np.isclose(I_dev, 1.0, rtol=1e-10)

        # Test scaling with Ie
        I_scaled, _ = sersic_profile(1.0, Ie=2.0)
        assert np.isclose(I_scaled, 2.0, rtol=1e-10)

    def test_sersic_profile_boundary_behavior(self) -> None:
        """Test behavior at boundaries."""
        # At R=0, profile should be maximum
        I_center, _ = sersic_profile(0.001)  # Very close to center
        I_outer, _ = sersic_profile(10.0)  # Far out

        # Center should be brighter than outer regions
        assert I_center > I_outer

    def test_get_luminosity_sersic_basic(self) -> None:
        """Test basic functionality of luminosity density calculation."""
        # Test that function runs without error
        j = get_luminosity_sersic(1.0)

        # Should return a float
        assert isinstance(j, float)

        # For reasonable parameters, should be positive
        assert j > 0

    def test_get_luminosity_sersic_scaling(self) -> None:
        """Test scaling behavior of luminosity density."""
        j1 = get_luminosity_sersic(1.0, Ie=1.0)
        j2 = get_luminosity_sersic(1.0, Ie=2.0)

        # Should scale linearly with Ie
        assert np.isclose(j2 / j1, 2.0, rtol=1e-2)

    @patch("matplotlib.pyplot.subplots")
    def test_plot_sersic_profile_creates_figure(self, mock_subplots: MagicMock) -> None:
        """Test that plot_sersic_profile creates a matplotlib figure."""

        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_subplots.return_value = (mock_fig, mock_ax)

        output_path = Path("docs/assets/generated/sersic_profile.png")
        plot_sersic_profile(output_path)

        # Verify figure creation
        mock_subplots.assert_called_once()
        mock_fig.subplots_adjust.assert_called_once()

        # Verify axes configuration
        mock_ax.set_xlabel.assert_called_once()
        mock_ax.set_ylabel.assert_called_once()
        mock_ax.set_xlim.assert_called_once()

        # Verify plot elements
        assert mock_ax.plot.call_count >= 1  # At least one profile plotted

        # Verify figure saving
        mock_fig.savefig.assert_called_once_with(
            output_path,
            dpi=150,
            bbox_inches="tight",
            facecolor="white",
            edgecolor="none",
        )

    @patch("matplotlib.pyplot.subplots")
    def test_plot_sersic_profile_shows_figure(self, mock_subplots: MagicMock) -> None:
        """Test that plot_sersic_profile shows figure when no path is provided."""
        with patch("galactic_dynamics_bt.chapter01.sersic_profile.plt.show") as mock_show:
            mock_fig = MagicMock()
            mock_ax = MagicMock()
            mock_subplots.return_value = (mock_fig, mock_ax)

            plot_sersic_profile()

            # Verify show was called and savefig was not called
            mock_show.assert_called_once()
            mock_fig.savefig.assert_not_called()

    def test_sersic_profile_mathematical_consistency(self) -> None:
        """Test mathematical consistency of the profile."""
        # Test that profile is monotonically decreasing
        radii = np.array([0.5, 1.0, 2.0, 5.0])
        intensities = [sersic_profile(r)[0] for r in radii]

        # Should be decreasing
        for i in range(len(intensities) - 1):
            assert intensities[i] > intensities[i + 1]
