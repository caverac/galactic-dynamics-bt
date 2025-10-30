# Contributing to Galactic Dynamics Solutions

Thank you for your interest in contributing to this project! This guide will help you get started with contributing solutions, improvements, and documentation to the Galactic Dynamics problem set solutions.

## 🎯 Ways to Contribute

- **Add new problem solutions** from Binney & Tremaine
- **Improve existing solutions** with better algorithms or explanations
- **Fix bugs** in calculations or code
- **Enhance documentation** with clearer explanations
- **Add visualizations** and interactive plots
- **Improve code quality** and performance
- **Write tests** for numerical solutions

## 🚀 Getting Started

### 1. Set Up Development Environment

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/your-username/galactic-dynamics-bt.git
   cd galactic-dynamics-bt
   ```

3. **Install uv** (if not already installed):
   ```bash
   curl -LsSf https://astral.sh/uv/install.sh | sh
   ```

4. **Set up the development environment**:
   ```bash
   uv sync --all-groups
   ```

5. **Install pre-commit hooks**:
   ```bash
   uv run pre-commit install
   ```

### 2. Create a Feature Branch

```bash
git checkout -b feature/chapter-X-problem-Y
# or
git checkout -b fix/issue-description
# or
git checkout -b docs/improve-explanation
```

## 📝 Contributing Guidelines

### Code Style and Quality

We maintain high code quality standards:

#### Python Code Style
- **Formatter**: [Black](https://black.readthedocs.io/) (line length: 120)
- **Linter**: [Flake8](https://flake8.pycqa.org/)
- **Type Checker**: [MyPy](https://mypy.readthedocs.io/)


Run quality checks:
```bash
# Format code
uv run black src tests

# Check formatting
uv run black --check src tests

# Lint code
uv run flake8 src

# Type checking
uv run mypy src

# Sort imports
uv run isort src tests
```

#### Documentation Style
- **Markdown**: Use [Material for MkDocs](https://squidfunk.github.io/mkdocs-material/) syntax
- **Math**: Use LaTeX syntax with custom macros (e.g., `\rmarcsec` for arcseconds)
- **Figures**: Save as PNG in `docs/assets/` with descriptive names

### Adding New Problem Solutions

#### 1. Code Structure

Create new solutions following this structure:

```python
"""Problem X.Y: Brief description.

Detailed explanation of the physical problem and approach.
"""

import logging
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

logger = logging.getLogger(__name__)


def solve_problem_xy() -> tuple[np.ndarray, np.ndarray]:
    """Solve Problem X.Y analytically/numerically.
    
    Returns:
        tuple: (x_values, y_values) for plotting
    """
    # Implementation here
    return x_values, y_values


def plot_problem_xy() -> None:
    """Generate publication-quality plot for Problem X.Y."""
    # Configure matplotlib for Times New Roman font
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times New Roman', 'Times']
    plt.rcParams['mathtext.fontset'] = 'stix'
    
    fig, ax = plt.subplots(figsize=(6.5, 5))
    
    # Your plotting code here
    
    # Save figure
    fig.savefig(
        f"docs/assets/chapter{X}_problem{Y}.png",
        dpi=300,
        bbox_inches='tight',
        facecolor='white',
        edgecolor='none'
    )


if __name__ == "__main__":
    plot_problem_xy()
```

#### 2. Documentation

Add corresponding documentation in `docs/chapterXX.md`:

```markdown
## Problem X.Y

Brief description of the problem setup.

### Physical Setup

Detailed explanation with equations:

$$
\frac{d^2 r}{dt^2} = -\frac{GM}{r^2}
$$

### Solution

Step-by-step derivation...

### Results

![Problem X.Y Results](assets/chapterX_problemY.png)

*Figure X.Y: Description of the plot and key insights.*
```

#### 3. Tests

Write unit tests for your solution:

```python
"""Tests for Problem X.Y solution."""

import pytest
import numpy as np
from chapter0X.problem_xy import solve_problem_xy


class TestProblemXY:
    """Test cases for Problem X.Y."""
    
    def test_solution_properties(self):
        """Test that solution has expected properties."""
        x, y = solve_problem_xy()
        
        # Test array properties
        assert len(x) == len(y)
        assert np.all(np.isfinite(x))
        assert np.all(np.isfinite(y))
        
    def test_physical_constraints(self):
        """Test that solution satisfies physical constraints."""
        # Add specific physics tests
        pass
```

### Commit Message Convention

We use [Conventional Commits](https://www.conventionalcommits.org/):

```
type(scope): description

feat(chapter1): add solution for Problem 1.5
fix(frw): correct sign error in acceleration equation
docs(readme): update installation instructions
test(chapter1): add tests for Problem 1.1
```

**Types:**
- `feat`: New feature or solution
- `fix`: Bug fix
- `docs`: Documentation changes
- `test`: Adding or updating tests
- `refactor`: Code refactoring
- `style`: Code style changes
- `chore`: Maintenance tasks

## 🧪 Testing

### Running Tests

```bash
# Run all tests
uv run pytest

# Run tests with coverage
uv run pytest --cov=src

# Run specific test file
uv run pytest tests/unit/chapter01/test_frw_model.py

# Run with verbose output
uv run pytest -v
```

### Writing Tests

- Write unit tests for all numerical functions
- Test edge cases and boundary conditions
- Verify physical constraints are satisfied
- Use descriptive test names and docstrings
- Mock matplotlib for plotting function tests

### Test Categories

1. **Unit Tests**: Test individual functions and methods
2. **Integration Tests**: Test complete problem solutions
3. **Physics Tests**: Verify solutions satisfy known physical laws
4. **Regression Tests**: Ensure numerical results remain stable

## 📚 Documentation Guidelines

### Mathematical Content

- Use clear, step-by-step derivations
- Include physical intuition and interpretation
- Reference relevant sections in Binney & Tremaine
- Use consistent notation throughout

### Figures and Plots

- Generate high-quality, publication-ready figures
- Use Times New Roman font for consistency
- Include descriptive captions
- Save as PNG with 300 DPI minimum
- Use meaningful filenames

### Code Documentation

- Write clear docstrings for all functions
- Include parameter types and return values
- Add inline comments for complex algorithms
- Provide usage examples

## 🔄 Pull Request Process

### Before Submitting

1. **Update your branch**:
   ```bash
   git fetch upstream
   git rebase upstream/main
   ```

2. **Run all quality checks**:
   ```bash
   uv run black src tests
   uv run flake8 src
   uv run mypy src
   uv run pytest
   ```

3. **Update documentation**:
   ```bash
   uv run mkdocs build --strict
   ```

4. **Commit your changes**:
   ```bash
   git add .
   git commit -m "feat(chapter1): add solution for Problem 1.5"
   ```

### Submitting the PR

1. **Push to your fork**:
   ```bash
   git push origin feature/chapter-X-problem-Y
   ```

2. **Create Pull Request** on GitHub with:
   - Clear title following conventional commits
   - Detailed description of changes
   - Reference to related issues
   - Screenshots/plots if applicable

3. **PR Template**:
   ```markdown
   ## Description
   Brief description of changes

   ## Problem/Issue
   - Solves Problem X.Y from Chapter X
   - Fixes issue #123

   ## Changes Made
   - [ ] Added analytical solution
   - [ ] Implemented numerical method
   - [ ] Created visualization
   - [ ] Added tests
   - [ ] Updated documentation

   ## Testing
   - [ ] All tests pass
   - [ ] Added new tests for functionality
   - [ ] Manual testing completed

   ## Screenshots/Plots
   (If applicable)
   ```

### Review Process

1. **Automated Checks**: CI runs tests and quality checks
2. **Code Review**: Maintainers review code and documentation
3. **Feedback**: Address any requested changes
4. **Approval**: Once approved, maintainers will merge

## 🎨 Style Guidelines

### Python Code

- Follow PEP 8 with 120 character line limit
- Use type hints for all function parameters and return values
- Use descriptive variable names (prefer `orbital_radius` over `r`)
- Group imports: standard library, third-party, local
- Use logging instead of print statements

### Mathematical Notation

- Follow Binney & Tremaine notation when possible
- Define custom LaTeX macros for repeated units
- Use consistent symbol definitions
- Include units in variable names and comments

### Plotting Style

```python
# Standard plot configuration
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'Times']
plt.rcParams['mathtext.fontset'] = 'stix'

# Tick configuration
ax.tick_params(direction='in', which='both')
ax.tick_params(top=True, right=True, which='both')
```

## 🛟 Getting Help

- **Questions**: Open a [GitHub Discussion](https://github.com/caverac/galactic-dynamics-bt/discussions)
- **Bugs**: Create an [Issue](https://github.com/caverac/galactic-dynamics-bt/issues)
- **Email**: Contact Carlos Vera-Ciro at [caverac@gmail.com](mailto:caverac@gmail.com)

## 📄 License

By contributing to this project, you agree that your contributions will be licensed under the MIT License.

Thank you for contributing to advancing galactic dynamics education! 🌌