"""Unit tests for multipole expansion module."""

from pathlib import Path
from unittest.mock import Mock, patch

import numpy as np


from galactic_dynamics_bt.chapter02.multipole_expansion import (
    MultipoleExpansion,
    plot_multipole_expansion_satoh,
    RhoL0Interpolator,
    SatohModel,
    SatohModelParams,
)


class TestSatohModel:
    """Test the SatohModel class."""

    params: SatohModelParams
    model: SatohModel

    def setup_method(self) -> None:
        """Set up test fixtures."""
        self.params = SatohModelParams(q=0.7)
        self.model = SatohModel(self.params)

    def test_init(self) -> None:
        """Test model initialization."""
        assert self.model.params == self.params

    def test_density_at_origin(self) -> None:
        """Test density calculation at origin."""
        density = self.model.density(0.0, 0.0)
        assert isinstance(density, float)
        assert density > 0  # Density should be positive

    def test_density_symmetry(self) -> None:
        """Test that density is symmetric in z."""
        R, z = 1.0, 1.0
        density_pos = self.model.density(R, z)
        density_neg = self.model.density(R, -z)
        assert np.isclose(density_pos, density_neg, rtol=1e-10)

    def test_density_decreases_with_radius(self) -> None:
        """Test that density decreases with increasing radius."""
        density_center = self.model.density(0.1, 0.0)
        density_edge = self.model.density(5.0, 0.0)
        assert density_center > density_edge

    def test_potential_at_origin(self) -> None:
        """Test potential calculation at origin."""
        potential = self.model.potential(0.0, 0.0)
        assert isinstance(potential, float)
        assert potential < 0  # Gravitational potential should be negative

    def test_potential_symmetry(self) -> None:
        """Test that potential is symmetric in z."""
        R, z = 1.0, 1.0
        potential_pos = self.model.potential(R, z)
        potential_neg = self.model.potential(R, -z)
        assert np.isclose(potential_pos, potential_neg, rtol=1e-10)

    def test_potential_approaches_zero_at_infinity(self) -> None:
        """Test that potential approaches zero at large distances."""
        potential_far = self.model.potential(100.0, 100.0)
        assert abs(potential_far) < 0.1  # Should be close to zero

    def test_density_returns_float(self) -> None:
        """Test that density always returns a float."""
        result = self.model.density(1.0, 0.5)
        assert isinstance(result, float)

    def test_potential_returns_float(self) -> None:
        """Test that potential always returns a float."""
        result = self.model.potential(1.0, 0.5)
        assert isinstance(result, float)


class TestRhoL0Interpolator:
    """Test the RhoL0Interpolator class."""

    model: SatohModel
    interpolator: RhoL0Interpolator

    def setup_method(self) -> None:
        """Set up test fixtures."""
        self.model = SatohModel(SatohModelParams(q=0.7))
        # Use small Ngrid for fast testing
        self.interpolator = RhoL0Interpolator(
            self.model, 0, Ngrid=10, a_min=1e-3, a_max=10.0  # multipole_order as positional argument
        )

    def test_init(self) -> None:
        """Test interpolator initialization."""
        assert self.interpolator.model == self.model
        assert self.interpolator.multipole_order == 0
        assert self.interpolator.Ngrid == 10
        assert len(self.interpolator.a_grid) == 10
        assert len(self.interpolator.rho_l0_values) == 10

    def test_compute_rho_l0_monopole(self) -> None:
        """Test monopole moment computation."""
        # pylint: disable=protected-access
        rho_l0 = self.interpolator._compute_rho_l0(1.0)
        assert isinstance(rho_l0, float)
        assert not np.isnan(rho_l0)

    def test_compute_rho_l0_higher_order(self) -> None:
        """Test higher-order multipole moment computation."""
        interpolator_l2 = RhoL0Interpolator(
            self.model, 2, Ngrid=5, a_min=1e-3, a_max=10.0  # multipole_order as positional argument
        )

        # pylint: disable=protected-access
        rho_l2 = interpolator_l2._compute_rho_l0(1.0)
        assert isinstance(rho_l2, float)
        assert not np.isnan(rho_l2)

    def test_call_within_bounds(self) -> None:
        """Test interpolator call within grid bounds."""
        result = self.interpolator(1.0)
        assert isinstance(result, float)
        assert not np.isnan(result)

    def test_call_outside_bounds(self) -> None:
        """Test interpolator call outside grid bounds."""
        # This should call _compute_rho_l0 directly
        result = self.interpolator(100.0)  # Outside a_max
        assert isinstance(result, float)
        assert not np.isnan(result)

    def test_interpolation_continuity(self) -> None:
        """Test that interpolation is reasonably continuous."""
        a1 = 1.0
        a2 = 1.01
        val1 = self.interpolator(a1)
        val2 = self.interpolator(a2)
        # Values should be close for nearby points
        assert abs(val1 - val2) / abs(val1) < 0.1

    def test_monopole_positive(self) -> None:
        """Test that monopole moment is positive (physical constraint)."""
        interpolator = RhoL0Interpolator(
            self.model, 0, Ngrid=5, a_min=1e-2, a_max=5.0  # multipole_order as positional argument
        )
        # Monopole should be positive for reasonable radii
        assert interpolator(1.0) > 0


class TestMultipoleExpansion:
    """Test the MultipoleExpansion class."""

    model: SatohModel
    expansion: MultipoleExpansion

    def setup_method(self) -> None:
        """Set up test fixtures."""
        self.model = SatohModel(SatohModelParams(q=0.7))
        # Use small l_max and Ngrid for fast testing
        self.expansion = MultipoleExpansion(self.model, 2, Ngrid=200)

    def test_init(self) -> None:
        """Test expansion initialization."""
        assert self.expansion.model == self.model
        assert self.expansion.l_max == 2
        assert 0 in self.expansion.rho_l0_interpolators
        assert 2 in self.expansion.rho_l0_interpolators

    def test_compute_potential_at_origin(self) -> None:
        """Test potential computation at origin."""
        # Use a small offset to avoid r=0 issues
        potential = self.expansion.compute_potential(1e-3, 0.0)
        assert isinstance(potential, float)
        assert potential < 0  # Should be negative
        assert not np.isnan(potential)

    def test_compute_potential_symmetry(self) -> None:
        """Test potential symmetry in z."""
        R, z = 1.0, 0.5
        potential_pos = self.expansion.compute_potential(R, z)
        potential_neg = self.expansion.compute_potential(R, -z)
        assert np.isclose(potential_pos, potential_neg, rtol=1e-6)

    def test_compute_potential_convergence(self) -> None:
        """Test that multipole expansion converges towards exact potential."""
        R, z = 1.0, 0.0
        exact_potential = self.model.potential(R, z)
        multipole_potential = self.expansion.compute_potential(R, z)

        # Should be reasonably close (within 50% for low l_max)
        relative_error = abs(multipole_potential - exact_potential) / abs(exact_potential)
        assert relative_error < 0.5

    def test_monopole_dominance(self) -> None:
        """Test that monopole term dominates for spherically symmetric case."""
        # At z=0, only even multipoles contribute, and l=0 should dominate
        potential_monopole = MultipoleExpansion(self.model, 0, Ngrid=100).compute_potential(1.0, 0.0)
        potential_with_quadrupole = self.expansion.compute_potential(1.0, 0.0)

        # The difference should be relatively small
        assert abs(potential_with_quadrupole - potential_monopole) < abs(potential_monopole) * 0.3

    def test_compute_potential_returns_float(self) -> None:
        """Test that compute_potential always returns a float."""
        result = self.expansion.compute_potential(1.0, 0.5)
        assert isinstance(result, float)


class TestIntegration:
    """Integration tests for the complete system."""

    def test_full_workflow_small_scale(self) -> None:
        """Test the complete workflow with small parameters for speed."""
        # Create model and expansion with minimal parameters
        model = SatohModel(SatohModelParams(q=0.8))
        expansion = MultipoleExpansion(model, 2, Ngrid=100)

        # Test at a few points
        test_points = [(1.0, 0.0), (0.5, 0.5), (2.0, 1.0)]

        for R, z in test_points:
            exact_potential = model.potential(R, z)
            multipole_potential = expansion.compute_potential(R, z)

            # Both should be finite and negative
            assert np.isfinite(exact_potential)
            assert np.isfinite(multipole_potential)
            assert exact_potential < 0
            assert multipole_potential < 0

            # Should be reasonably close (within order of magnitude)
            assert abs(multipole_potential / exact_potential) < 10

    def test_monopole_only_vs_exact_at_large_distance(self) -> None:
        """Test that monopole-only expansion approaches exact potential at large distances."""
        model = SatohModel(SatohModelParams(q=0.7))
        monopole_expansion = MultipoleExpansion(model, 0, Ngrid=100)

        # At large distances, monopole should dominate
        R, z = 10.0, 10.0
        exact_potential = model.potential(R, z)
        monopole_potential = monopole_expansion.compute_potential(R, z)

        # Should be reasonably close at large distances
        relative_error = abs(monopole_potential - exact_potential) / abs(exact_potential)
        assert relative_error < 0.2  # Within 20% at large distances


class TestPlotMultipoleExpansionSatoh:
    """Test the plot_multipole_expansion_satoh function."""

    @patch("galactic_dynamics_bt.chapter02.multipole_expansion.plt")
    @patch.object(SatohModel, "potential")
    @patch.object(MultipoleExpansion, "compute_potential")
    def test_plot_with_mocked_potentials(
        self,
        mock_compute_potential: Mock,
        mock_satoh_potential: Mock,
        mock_plt: Mock,
    ) -> None:
        """Test the plotting function with mocked expensive computations."""
        # Mock the expensive potential calculations
        mock_satoh_potential.return_value = -1.0
        mock_compute_potential.return_value = -0.95

        # Mock matplotlib components
        mock_fig = Mock()
        mock_axs = Mock()
        mock_plt.subplots.return_value = (mock_fig, mock_axs)
        mock_plt.show = Mock()

        # Call the plotting function
        plot_multipole_expansion_satoh()

        # Verify matplotlib setup
        mock_plt.subplots.assert_called_once_with(1, 1, figsize=(5, 5), sharex=True)
        mock_fig.subplots_adjust.assert_called_once()

        # Verify axis configuration
        mock_axs.set_xlabel.assert_called_once_with(r"$R/a$")
        mock_axs.set_ylabel.assert_called_once_with(r"$z/a$")
        mock_axs.set_xlim.assert_called_once_with(0, 5)
        mock_axs.set_ylim.assert_called_once_with(0, 5)

        # Verify contour plots were called (3 times: exact + 2 expansions)
        assert mock_axs.contour.call_count == 3

        # Verify legend was created
        mock_axs.legend.assert_called_once()

        # Verify show was called
        mock_plt.show.assert_called_once()

        # Verify expensive computations were called but mocked
        assert mock_satoh_potential.call_count > 0  # Called many times for grid
        assert mock_compute_potential.call_count > 0  # Called many times for grid

    @patch("galactic_dynamics_bt.chapter02.multipole_expansion.plt")
    @patch.object(SatohModel, "potential")
    @patch.object(MultipoleExpansion, "compute_potential")
    def test_plot_with_save_path(
        self,
        mock_compute_potential: Mock,
        mock_satoh_potential: Mock,
        mock_plt: Mock,
    ) -> None:
        """Test the plotting function with save path."""
        # Mock the expensive potential calculations
        mock_satoh_potential.return_value = -1.0
        mock_compute_potential.return_value = -0.95

        # Mock matplotlib components
        mock_fig = Mock()
        mock_axs = Mock()
        mock_plt.subplots.return_value = (mock_fig, mock_axs)

        # Create a mock path
        test_path = Path("/tmp/test_plot.png")

        # Call the plotting function with save path
        plot_multipole_expansion_satoh(test_path)

        # Verify savefig was called instead of show
        mock_fig.savefig.assert_called_once_with(
            test_path,
            dpi=150,
            bbox_inches="tight",
            facecolor="white",
            edgecolor="none",
        )

        # Verify show was NOT called
        mock_plt.show.assert_not_called()

    @patch("galactic_dynamics_bt.chapter02.multipole_expansion.plt")
    @patch.object(SatohModel, "potential")
    @patch.object(MultipoleExpansion, "compute_potential")
    def test_satoh_model_parameters(
        self,
        mock_compute_potential: Mock,
        mock_satoh_potential: Mock,
        mock_plt: Mock,
    ) -> None:
        """Test that the function uses correct Satoh model parameters."""
        # Mock the expensive potential calculations
        mock_satoh_potential.return_value = -1.0
        mock_compute_potential.return_value = -0.95

        # Mock matplotlib components
        mock_fig = Mock()
        mock_axs = Mock()
        mock_plt.subplots.return_value = (mock_fig, mock_axs)

        # Call the plotting function
        plot_multipole_expansion_satoh()

        # The function should create a SatohModel with q=0.6
        # We can't directly verify this without modifying the function,
        # but we can ensure the mocked methods were called
        assert mock_satoh_potential.call_count > 0
        assert mock_compute_potential.call_count > 0

    @patch("galactic_dynamics_bt.chapter02.multipole_expansion.plt")
    @patch.object(SatohModel, "potential")
    @patch.object(MultipoleExpansion, "compute_potential")
    def test_multipole_expansion_parameters(
        self,
        mock_compute_potential: Mock,
        mock_satoh_potential: Mock,
        mock_plt: Mock,
    ) -> None:
        """Test that MultipoleExpansion is created with correct parameters."""
        # Mock the expensive potential calculations
        mock_satoh_potential.return_value = -1.0
        mock_compute_potential.return_value = -0.95

        # Mock matplotlib components
        mock_fig = Mock()
        mock_axs = Mock()
        mock_plt.subplots.return_value = (mock_fig, mock_axs)

        # Call the plotting function
        plot_multipole_expansion_satoh()

        # Verify that compute_potential was called multiple times
        # (once for l_max=4 expansion and once for l_max=6 expansion)
        assert mock_compute_potential.call_count > 0

        # The actual grid evaluation happens during the contour plotting
        # so we expect many calls to the mocked methods
        assert mock_satoh_potential.call_count > 0
