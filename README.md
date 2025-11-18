# Galactic Dynamics (BT 2002) Solutions

[![Documentation](https://img.shields.io/badge/docs-mkdocs-blue.svg)](https://caverac.github.io/galactic-dynamics-bt/)
[![Python](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Solutions to the problem sets from Binney & Tremaine's *Galactic Dynamics* (2nd Edition). This repository contains Python implementations of numerical solutions, leveraging libraries such as NumPy, SciPy, and Matplotlib for efficient computation and visualization.

## 📚 About

This project provides detailed solutions to problems from the seminal textbook *Galactic Dynamics* by James Binney and Scott Tremaine (2nd Edition, 2008). The solutions are implemented in Python and organized by chapter, featuring:

- **Analytical derivations** with step-by-step mathematical explanations
- **Numerical implementations** using modern Python scientific computing stack
- **Visualizations** with publication-quality plots
- **Interactive documentation** with mathematical formulas and figures (**[https://caverac.github.io/galactic-dynamics-bt/](https://caverac.github.io/galactic-dynamics-bt/)**)

## 🚀 Quick Start

### Prerequisites

- Python 3.12 or higher
- [uv](https://github.com/astral-sh/uv) package manager (recommended) or pip

### Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/caverac/galactic-dynamics-bt.git
   cd galactic-dynamics-bt
   ```

2. Install dependencies using uv:
   ```bash
   uv sync
   ```

   Or using pip:
   ```bash
   pip install -e .
   ```

### Running the Code

Generate figures and run solutions:
```bash
# Generate plots for Chapter 1
uv run python scripts/generate_assets.py

# Run all tests
uv run pytest

# Build documentation
uv run mkdocs serve

# Build documentation - live
uv run mkdocs serve --livereload
```

## 📖 Documentation

The complete documentation with mathematical derivations, figures, and code explanations is available at:

**[https://caverac.github.io/galactic-dynamics-bt/](https://caverac.github.io/galactic-dynamics-bt/)**

### Local Documentation

To build and serve the documentation locally:

```bash
uv run mkdocs serve
```

Then open [http://localhost:8000](http://localhost:8000) in your browser.

## 📁 Project Structure

```
galactic-dynamics-bt/
├── docs/                          # Documentation source
│   ├── assets/                    # Generated figures and images
│   ├── chapterX.md                # Chapter X solutions
│   └── index.md                   # Homepage
├── galactic-dynamics-bt/          # Source code
│   └── chapterX/                  # Chapter X implementations
│       └── y.py                   # Implementation for problem y
├── tests/                         # Test suite
│   └── unit/                      # Unit tests
├── scripts/                       # Utility scripts
├── main.py                        # MkDocs macro definitions
├── pyproject.toml                 # Project configuration
└── mkdocs.yml                     # Documentation configuration
```

## 🔬 Implemented Solutions

### Chapter 1: Introduction

- **Problem 1.1**: Perihelion precession in constant density background
- **Problem 1.11**: FRW cosmological models and universe boundaries
- *More solutions coming soon...*

## 🧪 Development

### Setting up Development Environment

1. Install development dependencies:
   ```bash
   uv sync --group dev
   ```

2. Install pre-commit hooks:
   ```bash
   uv run pre-commit install
   ```

3. Run code quality checks:
   ```bash
   # Format code
   uv run black src tests

   # Lint code
   uv run flake8 src

   # Type checking
   uv run mypy src

   # Run tests
   uv run pytest
   ```

### Generating Releases

This project uses [Commitizen](https://commitizen-tools.github.io/commitizen/) for semantic versioning:

```bash
# Bump version and create changelog
uv run cz bump

# Push changes and tags
git push --follow-tags
```

## 🤝 Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for detailed guidelines on:

- Setting up the development environment
- Code style and quality standards
- Submitting pull requests
- Adding new problem solutions

## 📄 License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

## 🎯 Roadmap

- [x] Solutions for Chapter 1: Introduction
- [ ] Solutions for Chapter 2: Potential Theory
- [ ] Solutions for Chapter 3: The Orbits of Stars
- [ ] Solutions for Chapter 4: Equilibria of Collisionless Systems
- [ ] Solutions for Chapter 5: Stability of Collisionless Systems
- [ ] Solutions for Chapter 6: Disk Dynamics and Spiral Structure
- [ ] Solutions for Chapter 7: Kinetic Theory
- [ ] Solutions for Chapter 8: Collisions and Encounters of Stellar Systems
- [ ] Solutions for Chapter 9: Galaxy Formation
- [ ] Interactive Jupyter notebooks
- [ ] 3D visualizations for orbital mechanics
- [ ] Animation support for time-dependent solutions

## 📧 Contact

**Carlos Vera-Ciro** - [caverac@gmail.com](mailto:caverac@gmail.com)

Project Link: [https://github.com/caverac/galactic-dynamics-bt](https://github.com/caverac/galactic-dynamics-bt)

## 🙏 Acknowledgments

- James Binney and Scott Tremaine for their excellent textbook *Galactic Dynamics*
- The Python scientific computing community for the amazing tools
- Contributors and users of this repository

## 📚 References

Binney, J., & Tremaine, S. (2008). *Galactic Dynamics* (2nd ed.). Princeton University Press.
