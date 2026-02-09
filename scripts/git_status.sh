#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate bio-tools 2>/dev/null || true

cd /mnt/g/repos/bio
echo "=== git remote ==="
git remote -v

echo ""
echo "=== git status ==="
git status

echo ""
echo "=== git log (last 3) ==="
git log --oneline -3

echo ""
echo "=== untracked/new files ==="
git status --short
