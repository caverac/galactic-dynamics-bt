"""Oppenheimer-Volkoff equations for spherical equilibrium."""

import logging

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

logger = logging.getLogger(__name__)


def plot_frw_models() -> None:
    """Plot the FRW model parameter space."""

    fig: Figure
    axs: Axes
    fig, axs = plt.subplots(  # type: ignore[misc]
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
    y = 1 - x
    axs.plot(x, y, "-k", lw=1.5)
    axs.text(0.5, 0.4, "Open", ha="center", va="center", rotation=-40)
    axs.text(0.6, 0.5, "Closed", ha="center", va="center", rotation=-40)

    # accelerating-decelerating model
    x = np.linspace(0, 2.5, 100)
    y = 0.5 * x
    axs.plot(x, y, "--k", lw=1.5)

    fig.savefig(
        "docs/assets/frw_model.png",
        dpi=150,
        bbox_inches="tight",
        facecolor="white",
        edgecolor="none",
    )
