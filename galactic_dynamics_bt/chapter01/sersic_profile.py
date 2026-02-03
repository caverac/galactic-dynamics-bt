"""Sersic profile implementation and plotting for galaxy surface brightness distributions."""

from pathlib import Path

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import numpy as np
from scipy.integrate import quad

from galactic_dynamics_bt.utils.assets import register_asset, save_figure_if_changed


def sersic_profile(R: float, /, *, m: float = 4, Ie: float = 1, Re: float = 1) -> tuple[float, float]:
    """
    Calculate the surface brightness I(R) and its radial derivative dI(R)/dR for a Sersic profile.

    The Sersic profile is used to describe the surface brightness distribution of galaxies, and it
    follows the form:

        I(R) = Ie * exp(-b * ((R/Re)^(1/m) - 1))

    where b ~ 2*m - 0.324 is chosen so that Re contains half the total light, and m is the
    Sersic index. The classic de Vaucouleurs profile corresponds to m = 4.

    Parameters
    ----------
    R : float
        The projected radius at which to evaluate the profile, in the same units as Re.
    m : float, optional
        The Sersic index parameter. Default is 4 (classic de Vaucouleurs profile).
        - m = 1: Exponential profile (disk galaxies)
        - m = 4: de Vaucouleurs profile (elliptical galaxies)
        - m > 4: More concentrated profiles
    Ie : float, optional
        The surface brightness at the effective radius Re. Default is 1.
    Re : float, optional
        The effective radius (half-light radius) where the enclosed luminosity equals
        half the total luminosity. Default is 1.

    Returns
    -------
    I_R : float
        The surface brightness at radius R, in the same units as Ie.
    dI_dR : float
        The radial derivative of the surface brightness at radius R.


    Examples
    --------
    >>> sersic_profile(2.0)  # doctest: +ELLIPSIS
    (0.2340..., -0.2670...)

    """
    b = 2 * m - 0.324
    I_R = Ie * np.exp(-b * ((R / Re) ** (1 / m) - 1))

    dI_dR = -(b / (Re * m)) * (R / Re) ** (1 / m - 1) * I_R

    return float(I_R), float(dI_dR)


def get_luminosity_sersic(r: float, /, *, m: float = 4, Ie: float = 1, Re: float = 1) -> float:
    """
    Calculate the luminosity density j(r) for a Sersic profile.

    This function uses the Abel integral to deproject the surface brightness profile I(R)
    to obtain the three-dimensional luminosity density j(r). The deprojection assumes
    spherical symmetry.

    The Abel integral transforms the observed surface brightness I(R) into the
    intrinsic luminosity density j(r) via:

        j(r) = -(1/pi) integral[r to inf] (dI/dR) / sqrt(R^2 - r^2) dR

    Parameters
    ----------
    r : float
        The three-dimensional radius at which to calculate the luminosity density,
        in the same units as Re.
    m : float, optional
        The Sersic index parameter. Default is 4 (classic de Vaucouleurs profile).
        - m = 1: Exponential profile (disk galaxies)
        - m = 4: de Vaucouleurs profile (elliptical galaxies)
        - m > 4: More concentrated profiles
    Ie : float, optional
        The surface brightness at the effective radius Re. Default is 1.
    Re : float, optional
        The effective radius (half-light radius) where the enclosed luminosity equals
        half the total luminosity. Default is 1.

    Returns
    -------
    j_r : float
        The luminosity density at radius r, in units consistent with Ie and Re.

    Notes
    -----
    The Abel integral deprojection assumes that the galaxy has spherical symmetry.
    For real galaxies, which are often flattened, this provides an approximation
    to the intrinsic luminosity distribution.

    The integration is performed numerically using scipy.integrate.quad from r to
    infinity. Care should be taken when r approaches Re as the integrand may become
    steep near the effective radius.

    Examples
    --------
    >>> get_luminosity_sersic(2.0, m=4, Ie=1, Re=1.0)  # doctest: +ELLIPSIS
    0.0644803...

    """

    def integrand(R: float) -> float:
        _, dI_dR = sersic_profile(R, m=m, Ie=Ie, Re=Re)
        return float(dI_dR / np.sqrt(R * R - r * r))

    q, _ = quad(integrand, r, np.inf)

    return float(-q / np.pi)


@register_asset("sersic_profile.png")
def plot_sersic_profile(path: Path | None = None) -> None:
    """
    Plot the luminosity density profiles for different Sersic indices.

    Creates a log-log plot showing the deprojected luminosity density j(r)
    as a function of radius for three different Sersic indices (m = 3, 4, 5).
    The plot illustrates how the concentration of the profile changes with
    the Sersic index parameter.

    Parameters
    ----------
    path : Path | None, optional
        Path to save the figure, by default None. If None, the plot is shown instead.

    Notes
    -----
    - m = 3: Less concentrated profile
    - m = 4: Classic de Vaucouleurs profile (elliptical galaxies)
    - m = 5: More concentrated profile

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

    axs.set_xlabel(r"$\log_{10}r/R_e$")
    axs.set_ylabel(r"$\log_{10}j(r)R_e/I_e$")
    axs.set_xlim(-2, 2)

    # configure ticks
    axs.xaxis.set_major_locator(MultipleLocator(1))
    axs.yaxis.set_major_locator(MultipleLocator(1))

    axs.xaxis.set_minor_locator(AutoMinorLocator())
    axs.yaxis.set_minor_locator(AutoMinorLocator())

    axs.tick_params(which="minor", length=3, color="gray", direction="in")
    axs.tick_params(which="major", length=6, direction="in")
    axs.tick_params(top=True, right=True, which="both")

    # m = 4
    x = np.logspace(-2, 2, 100)

    y = np.array([get_luminosity_sersic(r, m=3, Ie=1, Re=1) for r in x])
    axs.plot(np.log10(x), np.log10(y), "--k", lw=1.5, label="$m=3$")

    y = np.array([get_luminosity_sersic(r, m=4, Ie=1, Re=1) for r in x])
    axs.plot(np.log10(x), np.log10(y), "-k", lw=1.5, label="$m=4$ (de Vaucouleurs)")

    y = np.array([get_luminosity_sersic(r, m=5, Ie=1, Re=1) for r in x])
    axs.plot(np.log10(x), np.log10(y), "-.k", lw=1.5, label="$m=5$")

    axs.legend(loc="lower left", frameon=False)

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
