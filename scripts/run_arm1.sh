#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate bio-tools 2>/dev/null || true
cd /mnt/g/repos/bio/bionemo
python 05_sstr2_smallmol_screen.py 2>&1
