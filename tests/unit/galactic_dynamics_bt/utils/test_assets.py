"""Test suite for assets utility module."""

from pathlib import Path
import tempfile

import matplotlib.pyplot as plt

from galactic_dynamics_bt.utils.assets import (
    figures_are_equal,
    get_registered_assets,
    save_figure_if_changed,
)


class TestGetRegisteredAssets:
    """Tests for get_registered_assets function."""

    def test_returns_dict(self) -> None:
        """Test that get_registered_assets returns a dictionary."""
        assets = get_registered_assets()
        assert isinstance(assets, dict)

    def test_returns_copy(self) -> None:
        """Test that get_registered_assets returns a copy, not the original."""
        assets1 = get_registered_assets()
        assets2 = get_registered_assets()
        # Should be equal but not the same object
        assert assets1 == assets2
        assert assets1 is not assets2


class TestFiguresAreEqual:
    """Tests for figures_are_equal function."""

    def test_returns_false_when_path_does_not_exist(self) -> None:
        """Test that figures_are_equal returns False when file doesn't exist."""
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3])

        result = figures_are_equal(Path("/nonexistent/path.png"), fig)

        assert result is False
        plt.close(fig)

    def test_returns_true_for_identical_figures(self) -> None:
        """Test that figures_are_equal returns True for identical figures."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "test.png"

            # Create and save a figure
            fig, ax = plt.subplots()
            ax.plot([1, 2, 3], [1, 4, 9])
            fig.savefig(path, dpi=100)

            # Compare with identical figure
            result = figures_are_equal(path, fig, dpi=100)

            assert result is True
            plt.close(fig)

    def test_returns_false_for_different_figures(self) -> None:
        """Test that figures_are_equal returns False for different figures."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "test.png"

            # Create and save a figure
            fig1, ax1 = plt.subplots()
            ax1.plot([1, 2, 3], [1, 4, 9])
            fig1.savefig(path, dpi=100)
            plt.close(fig1)

            # Create a different figure
            fig2, ax2 = plt.subplots()
            ax2.plot([1, 2, 3], [9, 4, 1])  # Different data

            result = figures_are_equal(path, fig2, dpi=100)

            assert result is False
            plt.close(fig2)

    def test_returns_false_for_different_sizes(self) -> None:
        """Test that figures_are_equal returns False for different figure sizes."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "test.png"

            # Create and save a figure with one size
            fig1, ax1 = plt.subplots(figsize=(4, 4))
            ax1.plot([1, 2, 3])
            fig1.savefig(path, dpi=100)
            plt.close(fig1)

            # Create a figure with different size
            fig2, ax2 = plt.subplots(figsize=(6, 6))
            ax2.plot([1, 2, 3])

            result = figures_are_equal(path, fig2, dpi=100)

            assert result is False
            plt.close(fig2)


class TestSaveFigureIfChanged:
    """Tests for save_figure_if_changed function."""

    def test_saves_when_file_does_not_exist(self) -> None:
        """Test that figure is saved when file doesn't exist."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "test.png"

            fig, ax = plt.subplots()
            ax.plot([1, 2, 3])

            result = save_figure_if_changed(fig, path)

            assert result is True
            assert path.exists()
            plt.close(fig)

    def test_skips_when_figure_unchanged(self) -> None:
        """Test that figure is skipped when unchanged."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "test.png"

            # Create and save initial figure
            fig, ax = plt.subplots()
            ax.plot([1, 2, 3], [1, 4, 9])
            save_figure_if_changed(fig, path, dpi=100)

            # Get mtime
            mtime_before = path.stat().st_mtime

            # Try to save same figure again
            result = save_figure_if_changed(fig, path, dpi=100)

            assert result is False
            # File should not have been modified
            assert path.stat().st_mtime == mtime_before
            plt.close(fig)

    def test_saves_when_figure_changed(self) -> None:
        """Test that figure is saved when it changed."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "test.png"

            # Create and save initial figure
            fig1, ax1 = plt.subplots()
            ax1.plot([1, 2, 3], [1, 4, 9])
            save_figure_if_changed(fig1, path, dpi=100)
            plt.close(fig1)

            # Create different figure
            fig2, ax2 = plt.subplots()
            ax2.plot([1, 2, 3], [9, 4, 1])

            result = save_figure_if_changed(fig2, path, dpi=100)

            assert result is True
            plt.close(fig2)

    def test_creates_parent_directories(self) -> None:
        """Test that parent directories are created if needed."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "subdir" / "nested" / "test.png"

            fig, ax = plt.subplots()
            ax.plot([1, 2, 3])

            result = save_figure_if_changed(fig, path)

            assert result is True
            assert path.exists()
            plt.close(fig)
