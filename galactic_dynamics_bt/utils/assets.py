"""Asset generation registry and utilities.

This module provides:

1. A decorator to register plot functions as documentation assets
2. A utility to save figures only when their content actually changes

Example usage:

    from galactic_dynamics_bt.utils.assets import register_asset

    @register_asset("chapter01/my_plot.png")
    def plot_my_figure(path: Path | None = None) -> None:
        ...
"""

from collections.abc import Callable
import io
from pathlib import Path
from typing import Any, Protocol, TypeVar

from matplotlib import pyplot as plt
from matplotlib.figure import Figure
import numpy as np


class PlotFunc(Protocol):
    """Protocol for plot functions that take an optional path keyword argument."""

    def __call__(self, path: Path | None = None) -> None:
        """Generate a plot, optionally saving to path."""


F = TypeVar("F", bound=PlotFunc)

# Global registry: output_name -> (function, module_name)
_registry: dict[str, tuple[PlotFunc, str]] = {}


def register_asset(output_name: str) -> Callable[[F], F]:
    """Define decorator to register a plot function as a documentation asset.

    Parameters
    ----------
    output_name : str
        Output filename relative to assets/generated/ directory.
        Example: "frw_models.png"

    Returns
    -------
    Callable
        Decorator that registers the function and returns it unchanged.

    Example
    -------
    >>> @register_asset("my_plot.png")
    ... def plot_something(path: Path | None = None) -> None:
    ...     fig, ax = plt.subplots()
    ...     ax.plot([1, 2, 3])
    ...     if path:
    ...         fig.savefig(path)
    ...     else:
    ...         plt.show()
    """

    def decorator(func: F) -> F:
        _registry[output_name] = (func, func.__module__)
        return func

    return decorator


def get_registered_assets() -> dict[str, tuple[PlotFunc, str]]:
    """Get all registered asset functions.

    Returns
    -------
    dict
        Mapping from output filename to (function, module_name) tuple.
    """
    return _registry.copy()


def figures_are_equal(path: Path, fig: Figure, **savefig_kwargs: Any) -> bool:
    """Check if a figure matches an existing file's pixel content.

    Parameters
    ----------
    path : Path
        Path to the existing image file.
    fig : Figure
        Matplotlib figure to compare.
    **savefig_kwargs
        Arguments passed to fig.savefig().

    Returns
    -------
    bool
        True if the figures have identical pixel content.
    """
    if not path.exists():
        return False

    # Render new figure to buffer
    buf = io.BytesIO()
    fig.savefig(buf, format="png", **savefig_kwargs)
    buf.seek(0)

    # Read both images as arrays
    new_img = plt.imread(buf)
    existing_img = plt.imread(path)

    # Compare pixel data
    if new_img.shape != existing_img.shape:
        return False

    return bool(np.array_equal(new_img, existing_img))


def save_figure_if_changed(
    fig: Figure,
    path: Path,
    **savefig_kwargs: Any,
) -> bool:
    """Save a figure only if its pixel content differs from the existing file.

    This prevents unnecessary file changes that would bloat git history
    when the visual content hasn't actually changed.

    Parameters
    ----------
    fig : Figure
        Matplotlib figure to save.
    path : Path
        Output path for the figure.
    **savefig_kwargs
        Arguments passed to fig.savefig() (e.g., dpi, bbox_inches).

    Returns
    -------
    bool
        True if the file was written, False if skipped (unchanged).
    """
    if figures_are_equal(path, fig, **savefig_kwargs):
        return False

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, **savefig_kwargs)
    return True
