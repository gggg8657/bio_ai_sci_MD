#!/usr/bin/env bash
# Push to GitHub: bio_ai_sci_MD
# Run: bash scripts/github_push.sh

set -e
cd /mnt/g/repos/bio

REPO_NAME="bio_ai_sci_MD"
GITHUB_USER="gggg8657"

echo "=== GitHub Push Script ==="

# Check if remote exists
if git remote get-url origin &>/dev/null; then
    echo "[INFO] Remote 'origin' already exists"
    git remote -v
else
    echo "[INFO] Adding remote origin..."
    git remote add origin "https://github.com/${GITHUB_USER}/${REPO_NAME}.git"
    echo "[OK] Remote added: https://github.com/${GITHUB_USER}/${REPO_NAME}.git"
fi

# Get current branch
BRANCH=$(git branch --show-current)
echo "[INFO] Current branch: $BRANCH"

# Push
echo "[INFO] Pushing to origin/$BRANCH..."
git push -u origin "$BRANCH"

echo "=== Done ==="
echo "Repository: https://github.com/${GITHUB_USER}/${REPO_NAME}"
