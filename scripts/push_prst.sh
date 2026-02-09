#!/usr/bin/env bash
set -e
cd /mnt/g/repos/bio

# Amend to include do_git.sh removal
git add -A
git diff --cached --quiet || git commit -m "chore: remove temp scripts"

# Add PRST_N_FM remote
git remote remove prst 2>/dev/null || true
git remote add prst https://github.com/MD-Agent/PRST_N_FM.git
echo "[OK] Remote added: https://github.com/MD-Agent/PRST_N_FM.git"

# Push to main
echo "[INFO] Pushing to prst/main..."
git push -u prst main --force
echo "[OK] Push complete"
echo "URL: https://github.com/MD-Agent/PRST_N_FM"
