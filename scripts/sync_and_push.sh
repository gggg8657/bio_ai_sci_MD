#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate bio-tools 2>/dev/null || true

cd /mnt/g/repos/bio

echo "=== Stash unstaged changes ==="
git stash 2>&1

echo ""
echo "=== Pull rebase from prst/main ==="
git pull prst main --rebase 2>&1

echo ""
echo "=== Pop stash ==="
git stash pop 2>&1

echo ""
echo "=== Push to prst ==="
git push prst main 2>&1

echo ""
echo "=== Final status ==="
git log --oneline -5
echo ""
git status --short
