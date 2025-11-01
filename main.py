"""MkDocs macros for dynamic version from pyproject.toml."""

from dataclasses import dataclass
from pathlib import Path
import tomllib
from typing import Any, Callable, Dict, Protocol

from dacite import from_dict


class MacroEnvironment(Protocol):
    """Type protocol for mkdocs-macros environment."""

    variables: Dict[str, Any]

    def macro(self, func: Callable[..., Any]) -> Callable[..., Any]:
        """Register a macro function."""
        ...


@dataclass
class PyProjectProject:
    """Data class for project section in pyproject.toml."""

    version: str


@dataclass
class PyprojectData:
    """Data class for pyproject.toml content."""

    project: PyProjectProject


def define_env(env: MacroEnvironment) -> None:
    """Define template variables for the mkdocs-macros plugin."""

    def get_project_version() -> str:
        """Get version from pyproject.toml."""
        pyproject_path = Path("pyproject.toml")
        if pyproject_path.exists():
            with open(pyproject_path, "rb") as f:
                data = from_dict(PyprojectData, tomllib.load(f))
                return data.project.version
        return "unknown"

    # Set the version as a variable that can be used in templates
    env.variables["version"] = get_project_version()

    # Also add it as a macro function
    @env.macro
    def project_version() -> str:  # pyright: ignore[reportUnusedFunction]
        """Return the project version."""
        return get_project_version()
