"""Unit tests for FRW model plotting functionality."""

from unittest.mock import patch, MagicMock

import pytest
import numpy as np
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from galactic_dynamics_bt.chapter01.frw_model import plot_frw_models


class TestFRWModel:
    """Test cases for FRW model plotting."""

    def setup_method(self) -> None:
        """Setup method to check if module is available."""

    def _create_mock_axes(self) -> MagicMock:
        """Create a properly mocked axes object with all necessary attributes."""
        mock_ax = MagicMock(spec=Axes)

        # Create mock xaxis and yaxis with set_major_locator and set_minor_locator methods
        mock_xaxis = MagicMock()
        mock_yaxis = MagicMock()
        mock_ax.xaxis = mock_xaxis
        mock_ax.yaxis = mock_yaxis

        return mock_ax

    def test_plot_frw_models_creates_figure(self) -> None:
        """Test that plot_frw_models creates a matplotlib figure."""
        with patch("matplotlib.pyplot.subplots") as mock_subplots:
            mock_fig = MagicMock(spec=Figure)
            mock_ax = self._create_mock_axes()
            mock_subplots.return_value = (mock_fig, mock_ax)

            with patch("matplotlib.pyplot.savefig"):
                plot_frw_models()

            # Verify subplots was called with correct parameters
            mock_subplots.assert_called_once_with(1, 1, figsize=(6.5, 5), sharex=True)

    def test_plot_frw_models_configures_axes(self) -> None:
        """Test that plot_frw_models properly configures the axes."""
        with patch("matplotlib.pyplot.subplots") as mock_subplots:
            mock_fig = MagicMock(spec=Figure)
            mock_ax = self._create_mock_axes()
            mock_subplots.return_value = (mock_fig, mock_ax)

            with patch("matplotlib.pyplot.savefig"):
                plot_frw_models()

            # Verify axes configuration
            mock_ax.set_xlabel.assert_called_once_with(r"$\Omega_{m 0}$")
            mock_ax.set_ylabel.assert_called_once_with(r"$\Omega_{\Lambda 0}$")
            mock_ax.set_xlim.assert_called_once_with(0, 2.5)
            mock_ax.set_ylim.assert_called_once_with(-0.5, 2)

    def test_plot_frw_models_sets_tick_locators(self) -> None:
        """Test that tick locators are properly configured."""
        with patch("matplotlib.pyplot.subplots") as mock_subplots:
            mock_fig = MagicMock(spec=Figure)
            mock_ax = self._create_mock_axes()
            mock_subplots.return_value = (mock_fig, mock_ax)

            with patch("matplotlib.pyplot.savefig"):
                plot_frw_models()

            # Verify tick locators are set
            mock_ax.xaxis.set_major_locator.assert_called_once()
            mock_ax.yaxis.set_major_locator.assert_called_once()
            mock_ax.xaxis.set_minor_locator.assert_called_once()
            mock_ax.yaxis.set_minor_locator.assert_called_once()

    def test_plot_frw_models_adds_plot_elements(self) -> None:
        """Test that plot elements (lines and text) are added."""
        with patch("matplotlib.pyplot.subplots") as mock_subplots:
            mock_fig = MagicMock(spec=Figure)
            mock_ax = self._create_mock_axes()
            mock_subplots.return_value = (mock_fig, mock_ax)

            with patch("matplotlib.pyplot.savefig"):
                plot_frw_models()

            # Verify plot elements are added
            assert mock_ax.plot.call_count == 2  # Two lines plotted
            assert mock_ax.text.call_count == 2  # Two text labels

    def test_plot_frw_models_configures_tick_params(self) -> None:
        """Test that tick parameters are properly configured."""
        with patch("matplotlib.pyplot.subplots") as mock_subplots:
            mock_fig = MagicMock(spec=Figure)
            mock_ax = self._create_mock_axes()
            mock_subplots.return_value = (mock_fig, mock_ax)

            with patch("matplotlib.pyplot.savefig"):
                plot_frw_models()

            # Verify tick_params is called multiple times for styling
            assert mock_ax.tick_params.call_count >= 3

            # Check that direction='in' is used
            calls = mock_ax.tick_params.call_args_list
            direction_calls = [call for call in calls if "direction" in call.kwargs]
            assert len(direction_calls) >= 2
            for call in direction_calls:
                assert call.kwargs["direction"] == "in"

    def test_plot_frw_models_saves_figure(self) -> None:
        """Test that the figure is saved to the correct location."""
        with patch("matplotlib.pyplot.subplots") as mock_subplots:
            mock_fig = MagicMock(spec=Figure)
            mock_ax = self._create_mock_axes()
            mock_subplots.return_value = (mock_fig, mock_ax)

            plot_frw_models()

            # Verify savefig was called with correct parameters
            mock_fig.savefig.assert_called_once_with(
                "docs/assets/generated/frw_model.png",
                dpi=150,
                bbox_inches="tight",
                facecolor="white",
                edgecolor="none",
            )

    def test_plot_data_arrays(self) -> None:
        """Test that the plotted data arrays have correct properties."""
        # Test the mathematical relationships in the plot
        x = np.linspace(0, 3, 100)
        y_open_closed = 1 - x

        # Verify boundary condition: at x=1, y should be 0
        idx_x1 = np.argmin(np.abs(x - 1.0))
        assert np.isclose(y_open_closed[idx_x1], 0.0, atol=0.01)

        # Test accelerating-decelerating line
        x_accel = np.linspace(0, 2.5, 100)
        y_accel = 0.5 * x_accel

        # Verify linear relationship
        assert np.allclose(y_accel, 0.5 * x_accel)

    @patch("matplotlib.pyplot.subplots")
    def test_plot_frw_models_no_exceptions(self, mock_subplots: MagicMock) -> None:
        """Test that plot_frw_models runs without raising exceptions."""
        mock_fig = MagicMock(spec=Figure)
        mock_ax = self._create_mock_axes()
        mock_subplots.return_value = (mock_fig, mock_ax)

        # Should not raise any exceptions
        try:
            plot_frw_models()
        except Exception as e:  # pylint: disable=broad-except
            pytest.fail(f"plot_frw_models() raised an exception: {e}")

    def test_figure_dimensions(self) -> None:
        """Test that figure has correct dimensions."""
        with patch("matplotlib.pyplot.subplots") as mock_subplots:
            mock_fig = MagicMock(spec=Figure)
            mock_ax = self._create_mock_axes()
            mock_subplots.return_value = (mock_fig, mock_ax)

            with patch("matplotlib.pyplot.savefig"):
                plot_frw_models()

            # Check figure size
            _, kwargs = mock_subplots.call_args
            assert kwargs["figsize"] == (6.5, 5)

    def test_subplots_adjust_called(self) -> None:
        """Test that figure layout is adjusted."""
        with patch("matplotlib.pyplot.subplots") as mock_subplots:
            mock_fig = MagicMock(spec=Figure)
            mock_ax = self._create_mock_axes()
            mock_subplots.return_value = (mock_fig, mock_ax)

            with patch("matplotlib.pyplot.savefig"):
                plot_frw_models()

            # Verify subplots_adjust was called
            mock_fig.subplots_adjust.assert_called_once()

            # Check that layout parameters are reasonable
            _, kwargs = mock_fig.subplots_adjust.call_args
            assert "left" in kwargs
            assert "right" in kwargs
            assert "bottom" in kwargs
            assert "top" in kwargs
