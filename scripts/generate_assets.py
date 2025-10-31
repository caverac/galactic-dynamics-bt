"""Generate asset files for documentation and testing."""

import os

from galactic_dynamics_bt.chapter01.frw_model import plot_frw_models
from galactic_dynamics_bt.chapter01.sersic_profile import plot_sersic_profile

if __name__ == "__main__":

    os.makedirs("docs/assets/generated", exist_ok=True)

    plot_frw_models()
    plot_sersic_profile()
