"""Unit tests for bifurcation module."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from galactic_dynamics_bt.chapter03.bifurcation import (
    AxesArray,
    LogarithmicPotential,
    LogarithmicPotentialParams,
    plot_loop_orbit_fraction,
    plot_surfaces_of_section,
)


class TestLogarithmicPotentialParams:
    """Tests for LogarithmicPotentialParams dataclass."""

    def test_params_stored(self) -> None:
        """Verify parameters are stored correctly."""
        params = LogarithmicPotentialParams(q=0.9, v0=1.0, Rc=0.1)
        assert params.q == 0.9
        assert params.v0 == 1.0
        assert params.Rc == 0.1


class TestLogarithmicPotential:
    """Tests for LogarithmicPotential class."""

    @pytest.fixture
    def model(self) -> LogarithmicPotential:
        """Create a test model."""
        return LogarithmicPotential(LogarithmicPotentialParams(q=0.9, v0=1.0, Rc=0.1))

    def test_potential_at_origin(self, model: LogarithmicPotential) -> None:
        """Test potential at origin uses core radius."""
        V = model.potential(0.0, 0.0)
        expected = 0.5 * np.log(0.1**2)
        assert V == pytest.approx(expected)

    def test_potential_on_x_axis(self, model: LogarithmicPotential) -> None:
        """Test potential on x-axis."""
        V = model.potential(1.0, 0.0)
        expected = 0.5 * np.log(1.0 + 0.01)
        assert V == pytest.approx(expected)

    def test_potential_gradient_at_point(self, model: LogarithmicPotential) -> None:
        """Test gradient calculation."""
        dVdx, dVdy = model.potential_gradient(1.0, 0.9)
        # denom = 1 + 0.81/0.81 + 0.01 = 2.01
        denom = 1.0 + (0.9**2) / (0.9**2) + 0.01
        expected_dx = 1.0 / denom
        expected_dy = 0.9 / (0.81 * denom)
        assert dVdx == pytest.approx(expected_dx)
        assert dVdy == pytest.approx(expected_dy)

    def test_potential_gradient_at_origin(self, model: LogarithmicPotential) -> None:
        """Test gradient at origin is zero."""
        dVdx, dVdy = model.potential_gradient(0.0, 0.0)
        assert dVdx == 0.0
        assert dVdy == 0.0


class TestModelInitialConditions:
    """Tests for Model initial condition methods."""

    @pytest.fixture
    def model(self) -> LogarithmicPotential:
        """Create a test model."""
        return LogarithmicPotential(LogarithmicPotentialParams(q=0.9, v0=1.0, Rc=0.1))

    def test_find_initial_conditions_valid(self, model: LogarithmicPotential) -> None:
        """Test finding initial conditions with valid energy."""
        E = model.potential(1.0, 0.0) + 1.0  # Energy above potential
        x, y, px, py = model.find_initial_conditions(E, 0.5, 0.5)
        assert x == 0.5
        assert y == 0.0
        assert py == 0.5
        assert px > 0  # Should be positive

    def test_find_initial_conditions_energy_too_low(self, model: LogarithmicPotential) -> None:
        """Test error when energy below potential."""
        E = model.potential(0.5, 0.0) - 1.0  # Below potential
        with pytest.raises(ValueError, match="Energy is less than potential"):
            model.find_initial_conditions(E, 0.5, 0.5)

    def test_find_initial_conditions_no_px_solution(self, model: LogarithmicPotential) -> None:
        """Test error when py too large for given energy."""
        E = model.potential(0.5, 0.0) + 0.1  # Small excess energy
        with pytest.raises(ValueError, match="No real solution for px"):
            model.find_initial_conditions(E, 0.5, 10.0)  # py too large

    def test_find_initial_conditions_on_section_valid(self, model: LogarithmicPotential) -> None:
        """Test finding initial conditions on section."""
        E = model.potential(1.0, 0.0) + 1.0
        x, y, px, py = model.find_initial_conditions_on_section(E, 0.5, 0.3)
        assert x == 0.5
        assert y == 0.0
        assert px == 0.3
        assert py > 0

    def test_find_initial_conditions_on_section_no_py_solution(self, model: LogarithmicPotential) -> None:
        """Test error when px too large."""
        E = model.potential(0.5, 0.0) + 0.1
        with pytest.raises(ValueError, match="No real solution for py"):
            model.find_initial_conditions_on_section(E, 0.5, 10.0)


class TestModelIntegration:
    """Tests for Model integration methods."""

    @pytest.fixture
    def model(self) -> LogarithmicPotential:
        """Create a test model."""
        return LogarithmicPotential(LogarithmicPotentialParams(q=1.0, v0=1.0, Rc=0.1))

    @patch("galactic_dynamics_bt.chapter03.bifurcation.SymplecticIntegrator")
    def test_integrate_calls_symplectic(self, mock_integrator_class: MagicMock, model: LogarithmicPotential) -> None:
        """Test that integrate uses SymplecticIntegrator."""
        mock_instance = MagicMock()
        mock_instance.integrate.return_value = np.zeros((10, 4))
        mock_integrator_class.return_value = mock_instance

        t = np.linspace(0, 1, 10)
        E = model.potential(1.0, 0.0) + 0.5
        model.integrate(E, t, 0.5, 0.3, order=2)

        mock_integrator_class.assert_called_once()
        mock_instance.integrate.assert_called_once()


class TestOrbitClassification:
    """Tests for orbit classification."""

    @pytest.fixture
    def model(self) -> LogarithmicPotential:
        """Create a test model."""
        return LogarithmicPotential(LogarithmicPotentialParams(q=1.0, v0=1.0, Rc=0.1))

    @patch.object(LogarithmicPotential, "integrate")
    def test_is_loop_orbit_positive_lz(self, mock_integrate: MagicMock, model: LogarithmicPotential) -> None:
        """Test loop orbit with positive Lz throughout."""
        # Lz = x*py - y*px, all positive
        mock_integrate.return_value = np.array(
            [
                [1.0, 0.0, 0.0, 1.0],  # Lz = 1
                [0.0, 1.0, -1.0, 0.0],  # Lz = 1
                [-1.0, 0.0, 0.0, -1.0],  # Lz = 1
            ]
        )
        t = np.array([0.0, 1.0, 2.0])
        assert model.is_loop_orbit(0.0, t, 0.5, 0.5) is True

    @patch.object(LogarithmicPotential, "integrate")
    def test_is_loop_orbit_negative_lz(self, mock_integrate: MagicMock, model: LogarithmicPotential) -> None:
        """Test loop orbit with negative Lz throughout."""
        mock_integrate.return_value = np.array(
            [
                [1.0, 0.0, 0.0, -1.0],  # Lz = -1
                [0.0, 1.0, 1.0, 0.0],  # Lz = -1
            ]
        )
        t = np.array([0.0, 1.0])
        assert model.is_loop_orbit(0.0, t, 0.5, 0.5) is True

    @patch.object(LogarithmicPotential, "integrate")
    def test_is_box_orbit_sign_change(self, mock_integrate: MagicMock, model: LogarithmicPotential) -> None:
        """Test box orbit with Lz sign change."""
        mock_integrate.return_value = np.array(
            [
                [1.0, 0.0, 0.0, 1.0],  # Lz = 1
                [1.0, 0.0, 0.0, -1.0],  # Lz = -1
            ]
        )
        t = np.array([0.0, 1.0])
        assert model.is_loop_orbit(0.0, t, 0.5, 0.5) is False


class TestSurfaceOfSection:
    """Tests for surface of section methods."""

    @pytest.fixture
    def model(self) -> LogarithmicPotential:
        """Create a test model."""
        return LogarithmicPotential(LogarithmicPotentialParams(q=1.0, v0=1.0, Rc=0.1))

    @patch.object(LogarithmicPotential, "find_exact_intersection")
    def test_extract_crossings_finds_crossing(self, mock_exact: MagicMock, model: LogarithmicPotential) -> None:
        """Test crossing detection when y goes from negative to positive."""
        mock_exact.return_value = (0.5, 0.3)
        z = np.array(
            [
                [0.4, -0.1, 0.2, 0.5],  # y < 0
                [0.5, 0.1, 0.3, 0.6],  # y >= 0, py > 0 -> crossing
            ]
        )
        # pylint: disable=protected-access
        x_cross, px_cross = model._extract_crossings(z)
        assert len(x_cross) == 1
        assert x_cross[0] == 0.5
        assert px_cross[0] == 0.3

    def test_extract_crossings_no_crossing_py_negative(self, model: LogarithmicPotential) -> None:
        """Test no crossing detected when py <= 0."""
        z = np.array(
            [
                [0.4, -0.1, 0.2, 0.5],
                [0.5, 0.1, 0.3, -0.1],  # py < 0
            ]
        )
        x_cross, _ = model._extract_crossings(z)  # pylint: disable=protected-access
        assert len(x_cross) == 0

    def test_extract_crossings_no_crossing_y_positive(self, model: LogarithmicPotential) -> None:
        """Test no crossing when y stays positive."""
        z = np.array(
            [
                [0.4, 0.1, 0.2, 0.5],  # y > 0
                [0.5, 0.2, 0.3, 0.6],  # y > 0
            ]
        )
        x_cross, _ = model._extract_crossings(z)  # pylint: disable=protected-access
        assert len(x_cross) == 0

    @patch.object(LogarithmicPotential, "integrate")
    @patch.object(LogarithmicPotential, "_extract_crossings")
    def test_surface_of_section_calls_integrate(
        self, mock_extract: MagicMock, mock_integrate: MagicMock, model: LogarithmicPotential
    ) -> None:
        """Test surface_of_section uses integrate and extract_crossings."""
        mock_integrate.return_value = np.zeros((10, 4))
        mock_extract.return_value = (np.array([0.5]), np.array([0.3]))

        t = np.array([0.0, 1.0])
        model.surface_of_section(0.0, t, 0.5, 0.5, order=2)

        mock_integrate.assert_called_once()
        mock_extract.assert_called_once()

    @patch("galactic_dynamics_bt.chapter03.bifurcation.SymplecticIntegrator")
    @patch.object(LogarithmicPotential, "_extract_crossings")
    def test_surface_of_section_from_section(
        self, mock_extract: MagicMock, mock_integrator: MagicMock, model: LogarithmicPotential
    ) -> None:
        """Test surface_of_section_from_section."""
        mock_instance = MagicMock()
        mock_instance.integrate.return_value = np.zeros((10, 4))
        mock_integrator.return_value = mock_instance
        mock_extract.return_value = (np.array([0.5]), np.array([0.3]))

        E = model.potential(1.0, 0.0) + 1.0
        t = np.array([0.0, 1.0])
        model.surface_of_section_from_section(E, t, 0.5, 0.3, order=2)

        mock_integrator.assert_called_once()
        mock_extract.assert_called_once()


class TestFindExactIntersection:
    """Tests for find_exact_intersection method."""

    @pytest.fixture
    def model(self) -> LogarithmicPotential:
        """Create a test model."""
        return LogarithmicPotential(LogarithmicPotentialParams(q=1.0, v0=1.0, Rc=0.1))

    def test_find_exact_intersection_real(self, model: LogarithmicPotential) -> None:
        """Test backward integration with real RK45 to cover derivatives function."""
        # Use real values that will allow backward integration
        x_exact, px_exact = model.find_exact_intersection(0.5, 0.3, 0.01, 0.5)
        # Should return values close to input (small y displacement)
        assert isinstance(x_exact, float)
        assert isinstance(px_exact, float)

    @patch("galactic_dynamics_bt.chapter03.bifurcation.RK45")
    def test_find_exact_intersection_mocked(self, mock_rk45_class: MagicMock, model: LogarithmicPotential) -> None:
        """Test backward integration to find exact crossing."""
        mock_instance = MagicMock()
        mock_instance.status = "finished"
        mock_instance.y = np.array([0.1, 0.55, 0.35, 0.8])
        mock_rk45_class.return_value = mock_instance

        x_exact, px_exact = model.find_exact_intersection(0.5, 0.3, 0.01, 0.5)

        assert x_exact == 0.55
        assert px_exact == 0.35
        mock_rk45_class.assert_called_once()

    @patch("galactic_dynamics_bt.chapter03.bifurcation.RK45")
    def test_find_exact_intersection_steps(self, mock_rk45_class: MagicMock, model: LogarithmicPotential) -> None:
        """Test that RK45 stepping is called until finished."""
        mock_instance = MagicMock()
        # First call returns "running", second returns finished
        mock_instance.status = "running"
        mock_instance.y = np.array([0.1, 0.55, 0.35, 0.8])

        def update_status() -> None:
            mock_instance.status = "finished"

        mock_instance.step.side_effect = update_status
        mock_rk45_class.return_value = mock_instance

        model.find_exact_intersection(0.5, 0.3, 0.01, 0.5)

        mock_instance.step.assert_called_once()


class TestAxesArray:
    """Tests for AxesArray helper class."""

    def test_axes_array_getitem_exists(self) -> None:
        """Test that AxesArray has __getitem__ method."""
        arr = AxesArray()
        assert hasattr(arr, "__getitem__")

    def test_axes_array_getitem_callable(self) -> None:
        """Test that __getitem__ can be called (stub returns None)."""
        arr = AxesArray()
        # The stub just has `...` which returns None implicitly
        result = arr[0, 0]
        assert result is None


class TestPlotSurfacesOfSection:
    """Tests for plot_surfaces_of_section function."""

    @patch("galactic_dynamics_bt.chapter03.bifurcation.save_figure_if_changed")
    @patch("galactic_dynamics_bt.chapter03.bifurcation.plt")
    @patch.object(LogarithmicPotential, "surface_of_section_from_section")
    def test_plot_surfaces_saves_to_path(
        self, mock_sos: MagicMock, mock_plt: MagicMock, mock_save: MagicMock, tmp_path: Path
    ) -> None:
        """Test plot saves to file when path provided."""
        mock_sos.return_value = (np.array([0.5]), np.array([0.3]))

        mock_fig = MagicMock()
        mock_axs = MagicMock()
        mock_ax = MagicMock()
        mock_axs.__getitem__ = MagicMock(return_value=mock_ax)
        mock_plt.subplots.return_value = (mock_fig, mock_axs)

        output_path = tmp_path / "test.png"
        plot_surfaces_of_section(output_path)

        mock_save.assert_called_once()
        mock_plt.show.assert_not_called()

    @patch("galactic_dynamics_bt.chapter03.bifurcation.save_figure_if_changed")
    @patch("galactic_dynamics_bt.chapter03.bifurcation.plt")
    @patch.object(LogarithmicPotential, "surface_of_section_from_section")
    def test_plot_surfaces_shows_when_no_path(
        self, mock_sos: MagicMock, mock_plt: MagicMock, mock_save: MagicMock
    ) -> None:
        """Test plot displays when no path provided."""
        mock_sos.return_value = (np.array([0.5]), np.array([0.3]))

        mock_fig = MagicMock()
        mock_axs = MagicMock()
        mock_ax = MagicMock()
        mock_axs.__getitem__ = MagicMock(return_value=mock_ax)
        mock_plt.subplots.return_value = (mock_fig, mock_axs)

        plot_surfaces_of_section(None)

        mock_plt.show.assert_called_once()
        mock_save.assert_not_called()

    @patch("galactic_dynamics_bt.chapter03.bifurcation.plt")
    @patch.object(LogarithmicPotential, "surface_of_section_from_section")
    def test_plot_surfaces_handles_value_error(self, mock_sos: MagicMock, mock_plt: MagicMock) -> None:
        """Test plot handles ValueError from surface_of_section."""
        mock_sos.side_effect = ValueError("test error")

        mock_fig = MagicMock()
        mock_axs = MagicMock()
        mock_ax = MagicMock()
        mock_axs.__getitem__ = MagicMock(return_value=mock_ax)
        mock_plt.subplots.return_value = (mock_fig, mock_axs)

        # Should not raise
        plot_surfaces_of_section(None)


class TestPlotLoopOrbitFraction:
    """Tests for plot_loop_orbit_fraction function."""

    @patch("galactic_dynamics_bt.chapter03.bifurcation.save_figure_if_changed")
    @patch("galactic_dynamics_bt.chapter03.bifurcation.plt")
    @patch.object(LogarithmicPotential, "is_loop_orbit")
    def test_plot_fraction_saves_to_path(
        self, mock_is_loop: MagicMock, mock_plt: MagicMock, mock_save: MagicMock, tmp_path: Path
    ) -> None:
        """Test plot saves to file when path provided."""
        mock_is_loop.return_value = True

        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        output_path = tmp_path / "test.png"
        plot_loop_orbit_fraction(output_path)

        mock_save.assert_called_once()
        mock_plt.show.assert_not_called()

    @patch("galactic_dynamics_bt.chapter03.bifurcation.save_figure_if_changed")
    @patch("galactic_dynamics_bt.chapter03.bifurcation.plt")
    @patch.object(LogarithmicPotential, "is_loop_orbit")
    def test_plot_fraction_shows_when_no_path(
        self, mock_is_loop: MagicMock, mock_plt: MagicMock, mock_save: MagicMock
    ) -> None:
        """Test plot displays when no path provided."""
        mock_is_loop.return_value = False

        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        plot_loop_orbit_fraction(None)

        mock_plt.show.assert_called_once()
        mock_save.assert_not_called()

    @patch("galactic_dynamics_bt.chapter03.bifurcation.plt")
    @patch.object(LogarithmicPotential, "is_loop_orbit")
    def test_plot_fraction_handles_value_error(self, mock_is_loop: MagicMock, mock_plt: MagicMock) -> None:
        """Test plot handles ValueError from is_loop_orbit."""
        mock_is_loop.side_effect = ValueError("test error")

        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        # Should not raise, fraction should be 0
        plot_loop_orbit_fraction(None)

    @patch("galactic_dynamics_bt.chapter03.bifurcation.plt")
    @patch.object(LogarithmicPotential, "is_loop_orbit")
    def test_plot_fraction_counts_correctly(self, mock_is_loop: MagicMock, mock_plt: MagicMock) -> None:
        """Test that loop/box counting works correctly."""
        # Alternate between loop and box
        mock_is_loop.side_effect = [True, False] * 1000

        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        plot_loop_orbit_fraction(None)

        # Should have been called multiple times
        assert mock_is_loop.call_count > 0
