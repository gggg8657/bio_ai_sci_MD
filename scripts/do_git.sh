#!/usr/bin/env bash
set -e
cd /mnt/g/repos/bio

echo "=== Checkout/create main branch ==="
git checkout -B main

echo "=== Stage all files ==="
git add -A

echo "=== Status ==="
git status --short

echo "=== Commit ==="
git commit -m "feat: restructure repo for PRST_N_FM with experiments and data

- Reorganize: data/fold_test1, results/foldmason, experiments/
- Add experiment docs: CIF->PDB, FoldMason MSA, PyMOL viz, PyRosetta setup
- Add README.md with project overview, pipeline diagram, quick start
- Include FoldMason results (HTML, FASTA, Newick)
- Update scripts paths for new directory structure
- Clean up redundant setup-only files"

echo "=== Log ==="
git log --oneline -3

echo "=== Done ==="
