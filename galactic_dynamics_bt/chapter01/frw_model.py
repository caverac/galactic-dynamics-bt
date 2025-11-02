"""Oppenheimer-Volkoff equations for spherical equilibrium."""

import logging
from pathlib import Path

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import numpy as np

logger = logging.getLogger(__name__)


def plot_frw_models(path: Path | None = None) -> None:
    """
    Plot the FRW model parameter space.

    Parameters
    ----------
    path : Path | None, optional
        Path to save the figure, by default None. If None, the plot is shown instead.
    """
    fig: Figure
    axs: Axes
    fig, axs = plt.subplots(
        1,
        1,
        figsize=(6.5, 5),
        sharex=True,
    )

    fig.subplots_adjust(
        left=0.1,
        right=0.92,
        bottom=0.1,
        top=0.92,
        wspace=0.05,
        hspace=0.05,
    )

    axs.set_xlabel(r"$\Omega_{m 0}$")
    axs.set_ylabel(r"$\Omega_{\Lambda 0}$")
    axs.set_xlim(0, 2.5)
    axs.set_ylim(-0.5, 2)

    # configure ticks
    axs.xaxis.set_major_locator(MultipleLocator(1))
    axs.yaxis.set_major_locator(MultipleLocator(1))

    axs.xaxis.set_minor_locator(AutoMinorLocator())
    axs.yaxis.set_minor_locator(AutoMinorLocator())

    axs.tick_params(which="minor", length=3, color="gray", direction="in")
    axs.tick_params(which="major", length=6, direction="in")
    axs.tick_params(top=True, right=True, which="both")

    # open-closed model
    x = np.linspace(0, 3, 100)
    y = np.ones_like(x) - x
    axs.plot(x, y, "-k", lw=1.5)
    axs.text(0.5, 0.4, "open", ha="center", va="center", rotation=-45, transform_rotates_text=True)
    axs.text(0.6, 0.5, "closed", ha="center", va="center", rotation=-45, transform_rotates_text=True)

    # accelerating-decelerating model
    x = np.linspace(0, 2.5, 100)
    y = 0.5 * x
    axs.plot(x, y, "--k", lw=1.5)
    axs.text(1.2, 0.7, "accelerating", ha="center", va="center", rotation=30, transform_rotates_text=True)
    axs.text(1.2, 0.5, "decelerating", ha="center", va="center", rotation=30, transform_rotates_text=True)

    # recollapse-expands forever model
    a = np.logspace(0.1, 4, 100)
    x = 2 * a**3 / (1 - 3 * a**2 + 2 * a**3)
    y = 1 / (1 - 3 * a**2 + 2 * a**3)
    axs.plot(x, y, "k", lw=1.5)
    axs.text(2.1, -0.02, "recollapses", ha="center", va="center", rotation=8, transform_rotates_text=True)
    axs.text(2.1, 0.14, "expands forever", ha="center", va="center", rotation=8, transform_rotates_text=True)

    # recollapse-expands forever model
    a = np.logspace(-1, -0.1, 100)
    x = 2 * a**3 / (1 - 3 * a**2 + 2 * a**3)
    y = 1 / (1 - 3 * a**2 + 2 * a**3)
    axs.plot(x, y, "k", lw=1.5)
    axs.fill_between(x, y.max() * np.ones_like(x), y, hatch=r"\\\\", edgecolor="black", facecolor="none")
    axs.text(0.3, 1.5, "bounce", ha="center", va="center", rotation=60, transform_rotates_text=True)

    if path:
        fig.savefig(
            path,
            dpi=150,
            bbox_inches="tight",
            facecolor="white",
            edgecolor="none",
        )
    else:
        plt.show()
