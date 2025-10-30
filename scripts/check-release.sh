#!/bin/bash

# Script to test Commitizen bump locally
# This helps verify the configuration before pushing to main

echo "🔍 Checking Commitizen configuration..."

# Check if commitizen is available
if ! uv run cz version > /dev/null 2>&1; then
    echo "❌ Commitizen not found. Run 'uv sync --group dev' first."
    exit 1
fi

echo "✅ Commitizen is available"

# Show current version
CURRENT_VERSION=$(uv run cz version --project)
echo "📌 Current version: $CURRENT_VERSION"

# Check if there are any conventional commits since last tag
echo "🔍 Checking for conventional commits..."

LAST_TAG=$(git describe --tags --abbrev=0 2>/dev/null || echo "")

if [ -z "$LAST_TAG" ]; then
    echo "ℹ️  No previous tags found"
    COMMITS=$(git log --oneline --grep="^feat" --grep="^fix" --grep="^perf" --grep="^refactor" --grep="^BREAKING CHANGE" --extended-regexp)
else
    echo "🏷️  Last tag: $LAST_TAG"
    COMMITS=$(git log ${LAST_TAG}..HEAD --oneline --grep="^feat" --grep="^fix" --grep="^perf" --grep="^refactor" --grep="^BREAKING CHANGE" --extended-regexp)
fi

if [ -n "$COMMITS" ]; then
    echo "✅ Found conventional commits:"
    echo "$COMMITS"
    echo ""
    echo "🚀 Would bump version. Run 'uv run cz bump --dry-run' to see what would happen."
else
    echo "ℹ️  No conventional commits found. Version bump not needed."
fi

echo ""
echo "💡 Commit message examples that trigger version bumps:"
echo "  feat: add new feature"
echo "  fix: fix bug"
echo "  perf: improve performance"
echo "  refactor: refactor code"
echo "  feat!: breaking change (bumps major version)"