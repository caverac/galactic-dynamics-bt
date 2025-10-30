# GitHub Actions Workflow Dependencies

This document explains the workflow dependencies for automated releases and documentation deployment.

## Workflow Sequence

When code is pushed to the `main` branch:

1. **Release Workflow** (`.github/workflows/release.yml`)
   - Triggers first on push to main
   - Checks for conventional commits
   - Bumps version in `pyproject.toml`
   - Updates `CHANGELOG.md`
   - Creates git tags
   - Pushes changes back to main
   - Creates GitHub release

2. **Documentation Workflow** (`.github/workflows/docs.yml`)
   - Triggers on workflow_run completion of the Release workflow
   - Also triggers on direct pushes (for cases without version bumps)
   - Fetches the latest changes (including updated version)
   - Builds documentation with the correct version
   - Deploys to GitHub Pages

## Key Features

### Release Workflow
- **Infinite Loop Prevention**: Uses `if: "!startsWith(github.event.head_commit.message, 'bump:')"` to prevent triggering on its own commits
- **Conditional Execution**: Only runs if conventional commits are found
- **Smart Version Detection**: Analyzes commits since last tag to determine appropriate version bump

### Documentation Workflow
- **Dependency Management**: Uses `workflow_run` trigger to wait for release completion
- **Fallback Trigger**: Still triggers on direct pushes for non-release changes
- **Fresh Checkout**: Fetches latest changes to ensure updated version is included

## Conventional Commit Examples

```bash
# Patch version bump (0.1.0 → 0.1.1)
git commit -m "fix: resolve calculation error in FRW model"

# Minor version bump (0.1.0 → 0.2.0) 
git commit -m "feat: add stellar orbit visualization"

# Major version bump (0.1.0 → 1.0.0)
git commit -m "feat!: redesign API for better usability"
```

## Testing Locally

Use the provided script to test release readiness:

```bash
./scripts/check-release.sh
```

This will show:
- Current version
- Conventional commits since last tag
- Whether a version bump would occur