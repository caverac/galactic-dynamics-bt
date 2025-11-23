"""Generate asset files for documentation and testing."""

import os
from pathlib import Path

from galactic_dynamics_bt.chapter01.frw_model import plot_frw_models
from galactic_dynamics_bt.chapter01.sersic_profile import plot_sersic_profile
from galactic_dynamics_bt.chapter01.universe_age import plot_universe_age
from galactic_dynamics_bt.chapter02.multipole_expansion import plot_multipole_expansion_satoh
from galactic_dynamics_bt.chapter03.surface_of_section import plot_angular_momentum, plot_orbital_properties

if __name__ == "__main__":

    path = Path("docs", "assets", "generated")
    os.makedirs(path, exist_ok=True)

    plot_frw_models(path / "frw_models.png")
    plot_sersic_profile(path / "sersic_profile.png")
    plot_universe_age(path / "universe_age.png")

    plot_multipole_expansion_satoh(path / "multipole_expansion_satoh.png")

    plot_orbital_properties(path / "orbital_properties.png")
    plot_angular_momentum(path / "angular_momentum.png")
