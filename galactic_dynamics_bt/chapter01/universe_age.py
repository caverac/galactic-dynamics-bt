"""Universe age calculations for FRW cosmological models."""

from pathlib import Path

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import numpy as np
from scipy.integrate import quad

from galactic_dynamics_bt.utils.assets import register_asset, save_figure_if_changed


def find_universe_age(
    z: float,
    /,
    *,
    omega_m0: float = 0.3,
    omega_lambda0: float = 0.7,
    omega_gamma0: float = 0.0,
    H0: float = 70.0,
) -> float:
    """Calculate the age of the universe in a flat FRW model.

    Parameters
    ----------
    omega_m0 : float
        The present-day matter density parameter.
    omega_lambda0 : float
        The present-day cosmological constant density parameter.
    omega_gamma0 : float, optional
        The present-day radiation density parameter. Default is 0.0.

    Returns
    -------
    float
        The age of the universe in units of Hubble time (1/H0).
    """

    def integrand(a: float) -> float:
        e_z = np.sqrt(
            omega_lambda0 + omega_m0 / a**3 + omega_gamma0 / a**4 + (1 - omega_m0 - omega_lambda0 - omega_gamma0) / a**2
        )
        return float(1 / (a * e_z))

    age, _ = quad(integrand, 0, 1 / (1 + z))
    return float(age * 1.02e3 / H0)


@register_asset("universe_age.png")
def plot_universe_age(path: Path | None = None) -> None:
    """
    Plot the age of the universe as a function of redshift.

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

    axs.set_xlabel(r"$1 + z$")
    axs.set_ylabel(r"$t/{\rm Gyr}$")
    axs.set_xscale("log")
    axs.set_yscale("log")
    axs.set_xlim(1000, 1)
    # axs.set_ylim(-0.5, 2)

    # configure ticks
    axs.tick_params(which="minor", length=3, color="gray", direction="in")
    axs.tick_params(which="major", length=6, direction="in")
    axs.tick_params(top=True, right=True, which="both")

    # open-closed model
    x = np.logspace(0, 3, 100)

    h7 = 1.05
    y = [
        find_universe_age(
            zp1 - 1,
            omega_m0=0.237,
            omega_lambda0=0.763,
            omega_gamma0=8.84e-5 / h7**2,
            H0=70.0 * h7,
        )
        for zp1 in x
    ]
    axs.plot(x, y, "-k", lw=1.5, label=r"$\Omega_{m0}=0.237$, $\Omega_{\Lambda0}=0.763$, $h=0.735$")

    # planck 2018 cosmology
    h = 0.674
    y = [
        find_universe_age(
            zp1 - 1,
            omega_m0=0.315,
            omega_lambda0=0.685,
            omega_gamma0=8.24e-5 / h**2,
            H0=100 * h,
        )
        for zp1 in x
    ]
    axs.plot(x, y, "--k", lw=1.5, label=r"$\Omega_{m0}=0.315$, $\Omega_{\Lambda0}=0.685$, $h=0.674$")

    axs.legend(loc="lower right", frameon=False)

    if path:
        save_figure_if_changed(
            fig,
            path,
            dpi=150,
            bbox_inches="tight",
            facecolor="white",
            edgecolor="none",
        )
    else:
        plt.show()
