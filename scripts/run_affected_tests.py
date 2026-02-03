#!/usr/bin/env python3
"""Run tests only for modified files.

This script maps source files to their corresponding test files and runs
only the affected tests. Used by pre-commit to avoid running the full
test suite on every commit.

Examples:
    galactic_dynamics_bt/chapter02/exponential_disk.py
        -> tests/unit/galactic_dynamics_bt/chapter02/test_exponential_disk.py

    tests/unit/galactic_dynamics_bt/chapter02/test_exponential_disk.py
        -> tests/unit/galactic_dynamics_bt/chapter02/test_exponential_disk.py
"""

from pathlib import Path
import subprocess
import sys

PACKAGE_NAME = "galactic_dynamics_bt"


def get_test_file_for_source(source_file: Path) -> Path | None:
    """Map a source file to its corresponding test file.

    Parameters
    ----------
    source_file : Path
        Path to the source or test file.

    Returns
    -------
    Path | None
        Path to the corresponding test file, or None if not found.
    """
    # Handle test files directly
    if source_file.name.startswith("test_"):
        return source_file

    # Only process Python files
    if source_file.suffix != ".py":
        return None

    parts = source_file.parts
    if PACKAGE_NAME not in parts:
        return None

    # Find the index of the package
    idx = parts.index(PACKAGE_NAME)
    relative_parts = parts[idx:]

    # Build test path
    test_path = Path("tests/unit") / Path(*relative_parts)
    test_path = test_path.with_name(f"test_{test_path.name}")

    if test_path.exists():
        return test_path

    return None


def main() -> int:
    """Run tests for modified files.

    Returns
    -------
    int
        Exit code (0 for success, non-zero for failure).
    """
    if len(sys.argv) < 2:
        print("No files provided")
        return 0

    files = [Path(f) for f in sys.argv[1:]]
    test_files: set[Path] = set()

    for file in files:
        test_file = get_test_file_for_source(file)
        if test_file and test_file.exists():
            test_files.add(test_file)

    if not test_files:
        print("No test files found for modified files")
        return 0

    print(f"Running tests: {', '.join(str(f) for f in sorted(test_files))}")

    # Run pytest without coverage to avoid the fail-under threshold
    result = subprocess.run(
        [sys.executable, "-m", "pytest", *[str(f) for f in test_files], "-v", "--no-cov"],
        check=False,
    )

    return result.returncode


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
