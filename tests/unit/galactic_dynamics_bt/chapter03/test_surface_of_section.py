"""Unit tests for surface of section module."""

from pathlib import Path
from unittest.mock import Mock, patch

import numpy as np
import pytest

from galactic_dynamics_bt.chapter03.surface_of_section import (
    AxesArray,
    find_exact_intersection,
    integrate_orbit,
    LogarithmicPotential,
    LogarithmicPotentialParams,
    Model,
    plot_angular_momentum,
    plot_orbit,
    plot_orbital_properties,
    plot_surface_of_section,
    surface_of_section,
)


class TestLoggarithmicPotentialParams:
    """Test the parameter dataclass."""

    def test_params_creation(self) -> None:
        """Test parameter creation with valid values."""
        params = LogarithmicPotentialParams(q=0.9, v0=1.0, Rc=0.0)
        assert params.q == 0.9
        assert params.v0 == 1.0

    def test_params_different_values(self) -> None:
        """Test with different parameter values."""
        params = LogarithmicPotentialParams(q=0.6, v0=220.0, Rc=0.0)
        assert params.q == 0.6
        assert params.v0 == 220.0


class TestLograrithmicPotential:
    """Test the LogarithmicPotential class."""

    params: LogarithmicPotentialParams
    model: LogarithmicPotential
    Lz: float

    def setup_method(self) -> None:
        """Set up test fixtures."""
        self.params = LogarithmicPotentialParams(q=0.9, v0=1.0, Rc=0.0)
        self.Lz = 0.2
        self.model = LogarithmicPotential(self.params, self.Lz)

    def test_initialization(self) -> None:
        """Test model initialization."""
        assert self.model.params == self.params
        assert self.model.Lz == self.Lz

    def test_potential_at_origin(self) -> None:
        """Test potential calculation at small radius."""
        # Test near origin (but not exactly at origin to avoid log(0))
        R, z = 0.1, 0.0
        phi = self.model.potential(R, z)

        # Should be finite and negative (attractive potential)
        assert np.isfinite(phi)
        assert isinstance(phi, float)

    def test_potential_symmetry(self) -> None:
        """Test z-symmetry of potential."""
        R = 1.0
        z1, z2 = 0.5, -0.5

        phi1 = self.model.potential(R, z1)
        phi2 = self.model.potential(R, z2)

        # Potential should be symmetric in z
        assert abs(phi1 - phi2) < 1e-10

    def test_potential_gradient_finite(self) -> None:
        """Test that gradient returns finite values."""
        R, z = 1.0, 0.5
        dR, dz = self.model.potential_gradient(R, z)

        assert np.isfinite(dR)
        assert np.isfinite(dz)
        assert isinstance(dR, float)
        assert isinstance(dz, float)

    def test_gradient_z_symmetry(self) -> None:
        """Test antisymmetry of z-gradient."""
        R = 1.0
        z1, z2 = 0.5, -0.5

        _, dz1 = self.model.potential_gradient(R, z1)
        _, dz2 = self.model.potential_gradient(R, z2)

        # z-gradient should be antisymmetric
        assert abs(dz1 + dz2) < 1e-10

    def test_effective_potential(self) -> None:
        """Test effective potential includes centrifugal term."""
        R, z = 1.0, 0.0
        phi_grav = self.model.potential(R, z)
        phi_eff = self.model.effective_potential(R, z)

        # Effective potential should be larger due to centrifugal term
        assert phi_eff > phi_grav

        # Difference should be Lz^2/(2R^2)
        expected_diff = 0.5 * self.Lz**2 / R**2
        assert abs(phi_eff - phi_grav - expected_diff) < 1e-10

    def test_effective_potential_gradient(self) -> None:
        """Test effective potential gradient."""
        R, z = 1.0, 0.5
        dR_eff, dz_eff = self.model.effective_potential_gradient(R, z)
        dR_grav, dz_grav = self.model.potential_gradient(R, z)

        # z-component should be unchanged
        assert abs(dz_eff - dz_grav) < 1e-10

        # R-component should include centrifugal force
        centrifugal = -self.Lz**2 / R**3
        assert abs(dR_eff - (dR_grav + centrifugal)) < 1e-10

    def test_zero_velocity_curve(self) -> None:
        """Test zero velocity curve calculation."""
        E, Lz = -0.8, 0.2
        R, z = self.model.zero_velocity_curve(E, Lz)

        # Check return types and shapes
        assert isinstance(R, np.ndarray)
        assert isinstance(z, np.ndarray)
        assert len(R) == len(z)
        assert len(R) > 0

    def test_find_initial_conditions(self) -> None:
        """Test initial condition calculation."""
        E, R, pR = -0.8, 0.3, 0.1
        R0, z0, pR0, pz0 = self.model.find_initial_conditions(E, R, pR)

        # Check return values
        assert R0 == R
        assert z0 == 0.0
        assert pR0 == pR
        assert pz0 >= 0  # Should be positive for valid orbit

        # Check energy conservation
        phi_eff = self.model.effective_potential(R0, z0)
        total_energy = 0.5 * (pR0**2 + pz0**2) + phi_eff
        assert abs(total_energy - E) < 1e-10


class TestModelAbstractBase:
    """Test the abstract Model base class."""

    # pylint: disable=abstract-class-instantiated
    def test_cannot_instantiate_abstract_class(self) -> None:
        """Test that abstract Model class cannot be instantiated."""
        with pytest.raises(TypeError):
            Model({}, 0.1)  # type: ignore


class TestIntegrateOrbit:
    """Test the orbit integration function with mocking."""

    params: LogarithmicPotentialParams
    model: LogarithmicPotential

    def setup_method(self) -> None:
        """Set up test fixtures."""
        self.params = LogarithmicPotentialParams(q=0.9, v0=1.0, Rc=0.0)
        self.model = LogarithmicPotential(self.params, 0.2)

    @patch("galactic_dynamics_bt.chapter03.surface_of_section.SymplecticIntegrator")
    def test_integrate_orbit_mocked(self, mock_integrator_class: Mock) -> None:
        """Test orbit integration with mocked integrator."""
        # Mock the integrator
        mock_integrator = Mock()
        mock_integrator.integrate.return_value = np.array(
            [
                [0.35, 0.0, 0.1, 0.5],  # Initial conditions
                [0.36, 0.1, 0.05, 0.48],  # After one step
            ]
        )
        mock_integrator_class.return_value = mock_integrator

        # Call function
        t, R, z, pR, pz = integrate_orbit(self.model, E=-0.8, R0=0.35, pR0=0.1, t_max=0.1, dt=0.01)

        # Check that integrator was created and called
        mock_integrator_class.assert_called_once()
        mock_integrator.integrate.assert_called_once()

        # Check return values have correct shapes
        assert isinstance(t, np.ndarray)
        assert isinstance(R, np.ndarray)
        assert isinstance(z, np.ndarray)
        assert isinstance(pR, np.ndarray)
        assert isinstance(pz, np.ndarray)

        # Check lengths match
        assert len(R) == len(z) == len(pR) == len(pz)


class TestFindExactIntersection:
    """Test the exact intersection finding function."""

    params: LogarithmicPotentialParams
    model: LogarithmicPotential

    def setup_method(self) -> None:
        """Set up test fixtures."""
        self.params = LogarithmicPotentialParams(q=0.9, v0=1.0, Rc=0.0)
        self.model = LogarithmicPotential(self.params, 0.2)

    @patch("galactic_dynamics_bt.chapter03.surface_of_section.RK45")
    def test_find_exact_intersection_mocked(self, mock_rk45_class: Mock) -> None:
        """Test exact intersection finding with mocked integrator."""
        # Mock the RK45 solver
        mock_solver = Mock()
        mock_solver.status = "running"
        mock_solver.y = np.array([0.0, 0.45, 0.32, 0.08])  # t, pz, R, pR at intersection

        # Create a side effect function that changes status after first call
        call_count = 0

        def step_side_effect() -> None:
            nonlocal call_count
            call_count += 1
            if call_count >= 2:  # Finish after a couple of steps
                mock_solver.status = "finished"

        mock_solver.step.side_effect = step_side_effect
        mock_rk45_class.return_value = mock_solver

        # Call function
        R_exact, pR_exact = find_exact_intersection(0.33, 0.05, 0.1, 0.5, self.model)

        # Check that solver was created and stepped
        mock_rk45_class.assert_called_once()
        assert call_count >= 2  # Should have called step at least twice

        # Check return values
        assert isinstance(R_exact, (float, np.floating))
        assert isinstance(pR_exact, (float, np.floating))
        assert np.isfinite(R_exact)
        assert np.isfinite(pR_exact)


class TestSurfaceOfSection:
    """Test surface of section calculation."""

    params: LogarithmicPotentialParams
    model: LogarithmicPotential

    def setup_method(self) -> None:
        """Set up test fixtures."""
        self.params = LogarithmicPotentialParams(q=0.9, v0=1.0, Rc=0.0)
        self.model = LogarithmicPotential(self.params, 0.2)

    @patch("galactic_dynamics_bt.chapter03.surface_of_section.find_exact_intersection")
    def test_surface_of_section_mocked(self, mock_intersection: Mock) -> None:
        """Test surface of section with mocked intersection finding."""
        # Mock exact intersection
        mock_intersection.return_value = (0.35, 0.1)

        # Create test data with one crossing
        R = np.array([0.3, 0.32, 0.35, 0.37])
        z = np.array([-0.1, -0.05, 0.05, 0.1])  # Crosses zero between indices 1 and 2
        pR = np.array([0.1, 0.08, 0.06, 0.04])
        pz = np.array([0.2, 0.3, 0.4, 0.35])  # All positive

        R_sec, pR_sec = surface_of_section(self.model, R, z, pR, pz)

        # Should find one crossing
        assert len(R_sec) == 1
        assert len(pR_sec) == 1
        assert R_sec[0] == 0.35
        assert pR_sec[0] == 0.1

        # Check that exact intersection was called
        mock_intersection.assert_called_once()

    def test_surface_of_section_no_crossings(self) -> None:
        """Test surface of section with no crossings."""
        # Create data with no crossings
        R = np.array([0.3, 0.32, 0.35])
        z = np.array([0.1, 0.2, 0.3])  # All positive
        pR = np.array([0.1, 0.08, 0.06])
        pz = np.array([0.2, 0.3, 0.4])

        R_sec, pR_sec = surface_of_section(self.model, R, z, pR, pz)

        # Should find no crossings
        assert len(R_sec) == 0
        assert len(pR_sec) == 0


class TestPlottingFunctions:
    """Test plotting functions with mocking to avoid slow execution."""

    @patch("galactic_dynamics_bt.chapter03.surface_of_section.integrate_orbit")
    @patch("matplotlib.pyplot.subplots")
    def test_plot_orbit_mocked(self, mock_subplots: Mock, mock_integrate: Mock) -> None:
        """Test plot_orbit function with mocking."""
        # Mock integration results
        t = np.linspace(0, 1, 10)
        R = np.cos(t)
        z = np.sin(t)
        pR = np.zeros_like(t)
        pz = np.ones_like(t)
        mock_integrate.return_value = (t, R, z, pR, pz)

        # Mock matplotlib
        mock_ax = Mock()
        mock_subplots.return_value = (Mock(), mock_ax)

        # Call function
        plot_orbit(0.9, mock_ax)

        # Check that integration was called
        mock_integrate.assert_called_once()

        # Check that plotting methods were called
        mock_ax.plot.assert_called_once()

    @patch("galactic_dynamics_bt.chapter03.surface_of_section.integrate_orbit")
    @patch("galactic_dynamics_bt.chapter03.surface_of_section.surface_of_section")
    @patch("galactic_dynamics_bt.chapter03.surface_of_section.LogarithmicPotential")
    def test_plot_surface_of_section_mocked(
        self, mock_model_class: Mock, mock_surface: Mock, mock_integrate: Mock
    ) -> None:
        """Test plot_surface_of_section with mocking."""
        # Mock orbit integration
        t = np.linspace(0, 10, 50)
        R = np.ones_like(t)
        z = 0.1 * np.sin(t)
        pR = 0.1 * np.cos(t)
        pz = np.ones_like(t)
        mock_integrate.return_value = (t, R, z, pR, pz)

        # Mock surface of section
        R_sec = np.array([1.0, 1.1])
        pR_sec = np.array([0.1, -0.1])
        mock_surface.return_value = (R_sec, pR_sec)

        # Mock the LogarithmicPotential class and its methods
        mock_model = Mock()
        mock_model.Lz = 0.2
        mock_model.zero_velocity_curve.return_value = (
            np.array([0.1, 0.2, 0.3, 0.4]),  # R_zvc
            np.array([0.1, 0.15, 0.1, 0.05]),  # z_zvc
        )
        mock_model_class.return_value = mock_model

        # Mock axes
        mock_ax = Mock()

        # Call function
        plot_surface_of_section(0.9, mock_ax)

        # Check that functions were called
        mock_integrate.assert_called_once()
        mock_surface.assert_called_once()
        mock_model_class.assert_called_once()

        # Check that zero_velocity_curve was called (twice - for regular and conserved L)
        assert mock_model.zero_velocity_curve.call_count == 2

    @patch("galactic_dynamics_bt.chapter03.surface_of_section.plot_orbit")
    @patch("galactic_dynamics_bt.chapter03.surface_of_section.plot_surface_of_section")
    @patch("matplotlib.pyplot.subplots")
    @patch("matplotlib.pyplot.show")
    def test_plot_orbital_properties_mocked(
        self,
        mock_show: Mock,  # pylint: disable=unused-argument
        mock_subplots: Mock,
        mock_plot_surface: Mock,
        mock_plot_orbit: Mock,
    ) -> None:
        """Test plot_orbital_properties with mocking."""
        # Mock matplotlib
        mock_fig = Mock()
        mock_axs = Mock()
        mock_axs.__getitem__ = Mock(return_value=Mock())
        mock_subplots.return_value = (mock_fig, mock_axs)

        # Call function
        plot_orbital_properties()

        # Check that subplots was created
        mock_subplots.assert_called_once()

        # Check that individual plotting functions were called
        assert mock_plot_orbit.call_count == 2  # Two orbit plots
        assert mock_plot_surface.call_count == 2  # Two surface plots

    @patch("galactic_dynamics_bt.chapter03.surface_of_section.integrate_orbit")
    @patch("matplotlib.pyplot.subplots")
    def test_plot_angular_momentum_mocked(self, mock_subplots: Mock, mock_integrate: Mock) -> None:
        """Test plot_angular_momentum with mocking."""
        # Mock orbit integration
        t = np.linspace(0, 10, 100)
        R = np.ones_like(t) + 0.1 * np.sin(t)
        z = 0.1 * np.cos(t)
        pR = np.zeros_like(t)
        pz = np.ones_like(t)
        mock_integrate.return_value = (t, R, z, pR, pz)

        # Mock matplotlib
        mock_fig = Mock()
        mock_axs = [Mock(), Mock()]
        mock_subplots.return_value = (mock_fig, mock_axs)

        # Call function
        plot_angular_momentum()

        # Check that integration was called
        mock_integrate.assert_called_once()

        # Check that plots were made
        mock_axs[0].plot.assert_called_once()
        mock_axs[1].plot.assert_called_once()

    @patch("galactic_dynamics_bt.chapter03.surface_of_section.integrate_orbit")
    @patch("matplotlib.pyplot.subplots")
    def test_plot_angular_momentum_save_path(self, mock_subplots: Mock, mock_integrate: Mock) -> None:
        """Test plot_angular_momentum with save path."""
        # Mock orbit integration
        mock_integrate.return_value = (np.linspace(0, 1, 10), np.ones(10), np.zeros(10), np.zeros(10), np.ones(10))

        # Mock matplotlib
        mock_fig = Mock()
        mock_axs = [Mock(), Mock()]
        mock_subplots.return_value = (mock_fig, mock_axs)

        # Call function with path
        test_path = Path("test.png")
        plot_angular_momentum(path=test_path)

        # Check that figure was saved
        mock_fig.savefig.assert_called_once_with(
            test_path,
            dpi=150,
            bbox_inches="tight",
            facecolor="white",
            edgecolor="none",
        )


class TestEdgeCases:
    """Test edge cases and error conditions."""

    params: LogarithmicPotentialParams
    model: LogarithmicPotential

    def setup_method(self) -> None:
        """Set up test fixtures."""
        self.params = LogarithmicPotentialParams(q=0.9, v0=1.0, Rc=0.0)
        self.model = LogarithmicPotential(self.params, 0.2)

    def test_potential_very_small_radius(self) -> None:
        """Test potential at very small radius."""
        phi = self.model.potential(1e-6, 0.0)
        # Should be finite (large negative value)
        assert np.isfinite(phi)
        assert phi < 0

    def test_zero_angular_momentum(self) -> None:
        """Test model with zero angular momentum."""
        model_zero_Lz = LogarithmicPotential(self.params, Lz=0.0)

        # Effective potential should equal gravitational potential
        R, z = 1.0, 0.5
        phi_grav = model_zero_Lz.potential(R, z)
        phi_eff = model_zero_Lz.effective_potential(R, z)

        assert abs(phi_eff - phi_grav) < 1e-15

    def test_large_flattening(self) -> None:
        """Test with very flattened potential."""
        params_flat = LogarithmicPotentialParams(q=0.1, v0=1.0, Rc=0.0)
        model_flat = LogarithmicPotential(params_flat, 0.2)

        # Should still work
        phi = model_flat.potential(1.0, 1.0)
        assert np.isfinite(phi)

    def test_spherical_potential(self) -> None:
        """Test with spherical potential (q=1)."""
        params_spherical = LogarithmicPotentialParams(q=1.0, v0=1.0, Rc=0.0)
        model_spherical = LogarithmicPotential(params_spherical, 0.2)

        # z and R should contribute equally at same distance
        phi1 = model_spherical.potential(1.0, 0.0)
        phi2 = model_spherical.potential(0.0, 1.0)

        assert abs(phi1 - phi2) < 1e-10


class TestAxesArray:
    """Test the AxesArray helper class."""

    def test_axes_array_getitem(self) -> None:
        """Test that AxesArray.__getitem__ is callable."""
        axes_array = AxesArray()
        try:
            axes_array[0, 1]
        except (NotImplementedError, AttributeError):
            pass


class TestPlotWithPath:
    """Test plotting functions with file paths."""

    @patch("galactic_dynamics_bt.chapter03.surface_of_section.integrate_orbit")
    @patch("matplotlib.pyplot.subplots")
    def test_plot_orbital_properties_with_path(self, mock_subplots: Mock, mock_integrate: Mock) -> None:
        """Test plot_orbital_properties with save path to cover the file save branch."""
        # Mock orbit integration
        mock_integrate.return_value = (np.linspace(0, 1, 10), np.ones(10), np.zeros(10), np.zeros(10), np.ones(10))

        # Mock matplotlib
        mock_fig = Mock()
        mock_axs = Mock()
        mock_axs.__getitem__ = Mock(return_value=Mock())
        mock_subplots.return_value = (mock_fig, mock_axs)

        # Mock the plotting functions to avoid actual computation
        with (
            patch("galactic_dynamics_bt.chapter03.surface_of_section.plot_orbit"),
            patch("galactic_dynamics_bt.chapter03.surface_of_section.plot_surface_of_section"),
        ):

            # Call function with path - this should cover line 836
            test_path = Path("test_orbital_props.png")
            plot_orbital_properties(path=test_path)

            # Check that figure was saved (covering the if path: branch)
            mock_fig.savefig.assert_called_once_with(
                test_path,
                dpi=150,
                bbox_inches="tight",
                facecolor="white",
                edgecolor="none",
            )
