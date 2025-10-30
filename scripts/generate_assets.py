"""Generate asset files for documentation and testing."""

import os
from galactic_dynamics_bt.chapter01.frw_model import plot_frw_models

if __name__ == "__main__":

    os.makedirs("docs/assets/generated", exist_ok=True)

    plot_frw_models()
