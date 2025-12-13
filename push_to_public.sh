#!/bin/bash
# push_to_public.sh
# Push to ExaPsi/QuasiX public repo, excluding private directories

set -e

# Directories and files to EXCLUDE from public repo
EXCLUDE_PATTERNS=(
    ".claude"
    "CLAUDE.md"
    "analysis"
    "output"
    "examples"
    "scripts"
    "tests"
    "verification_outputs"
    "docs"
    "backup"
)

# Build exclude arguments for git archive
EXCLUDE_ARGS=""
for pattern in "${EXCLUDE_PATTERNS[@]}"; do
    EXCLUDE_ARGS="$EXCLUDE_ARGS --exclude=$pattern"
done

echo "=== Pushing to ExaPsi/QuasiX (public) ==="
echo "Excluding: ${EXCLUDE_PATTERNS[*]}"
echo ""

# Create a temporary branch with filtered content
TEMP_BRANCH="public-filtered-$(date +%s)"
CURRENT_BRANCH=$(git branch --show-current)

echo "1. Creating filtered branch: $TEMP_BRANCH"

# Create orphan branch
git checkout --orphan $TEMP_BRANCH

# Remove all files from index
git rm -rf --cached . > /dev/null 2>&1 || true

# Add all files except excluded ones
git add -A

# Remove excluded directories from staging
for pattern in "${EXCLUDE_PATTERNS[@]}"; do
    git reset HEAD -- "$pattern" 2>/dev/null || true
    git rm -rf --cached "$pattern" 2>/dev/null || true
done

echo "2. Committing filtered content"
git commit -m "Public release - $(date +%Y-%m-%d)"

echo "3. Pushing to public remote"
git push public $TEMP_BRANCH:main --force

echo "4. Cleaning up"
git checkout $CURRENT_BRANCH
git branch -D $TEMP_BRANCH

echo ""
echo "=== Done! Pushed to ExaPsi/QuasiX (excluding private content) ==="
