#!/usr/bin/env python
# pylint: disable=unused-import
# flake8: noqa=F401
"""Generate documentation assets from solution plots.

This script discovers all registered plot functions and generates their
corresponding figure assets for the documentation.

Usage:
    uv run python scripts/generate_assets.py
"""

import re
from pathlib import Path

# Import chapter modules to trigger decorator registration
import galactic_dynamics_bt.chapter01.frw_model
import galactic_dynamics_bt.chapter01.sersic_profile
import galactic_dynamics_bt.chapter01.universe_age
import galactic_dynamics_bt.chapter02.multipole_expansion
import galactic_dynamics_bt.chapter03.surface_of_section
import galactic_dynamics_bt.chapter03.bifurcation
from galactic_dynamics_bt.utils.assets import get_registered_assets

ASSETS_DIR = Path(__file__).parent.parent / "docs" / "assets" / "generated"


def extract_chapter_number(module_name: str) -> tuple[int, int]:
    """Extract chapter and suborder from module name for sorting.

    Parameters
    ----------
    module_name : str
        Module name like 'galactic_dynamics_bt.chapter01.frw_model'

    Returns
    -------
    tuple[int, int]
        (chapter_number, suborder) for sorting. Suborder is based on
        the module name to ensure consistent ordering within chapters.
    """
    match = re.search(r"chapter(\d+)", module_name)
    chapter = int(match.group(1)) if match else 99
    # Use hash of module name for consistent sub-ordering within chapter
    suborder = hash(module_name) % 10000
    return (chapter, suborder)


def main() -> None:
    """Generate all documentation assets."""
    ASSETS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Assets directory: {ASSETS_DIR}")

    assets = get_registered_assets()
    print(f"Found {len(assets)} registered assets\n")

    written = 0
    skipped = 0

    # Sort by chapter number extracted from module path
    sorted_assets = sorted(
        assets.items(),
        key=lambda item: (extract_chapter_number(item[1][1]), item[0]),
    )

    current_chapter = None
    for output_name, (plot_func, module_name) in sorted_assets:
        # Print chapter header when chapter changes
        chapter_num = extract_chapter_number(module_name)[0]
        if chapter_num != current_chapter:
            current_chapter = chapter_num
            print(f"Chapter {current_chapter:02d}:")

        output_path = ASSETS_DIR / output_name

        # Check if file exists and get mtime before generation
        existed = output_path.exists()
        mtime_before = output_path.stat().st_mtime if existed else None

        plot_func(path=output_path)

        # Check if file was actually written
        if not existed or output_path.stat().st_mtime != mtime_before:
            print(f"  [WRITE] {output_name}")
            written += 1
        else:
            print(f"  [SKIP]  {output_name} (unchanged)")
            skipped += 1

    print(f"\nDone: {written} written, {skipped} unchanged")


if __name__ == "__main__":  # pragma: no cover
    main()
