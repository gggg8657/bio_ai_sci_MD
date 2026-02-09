#!/bin/bash
# SSTR2 Virtual Screening Pipeline - All 3 Arms
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate bio-tools 2>/dev/null || true

cd /mnt/g/repos/bio/bionemo

echo "========================================="
echo "  SSTR2 Virtual Screening Pipeline"
echo "========================================="
echo ""

echo "===== Arm 1: Small Molecule Screening ====="
python 05_sstr2_smallmol_screen.py 2>&1

echo ""
echo "===== Arm 2: FlexPepDock ====="
python 06_sstr2_flexpep_dock.py 2>&1

echo ""
echo "===== Arm 3: De Novo Binder Design ====="
python 07_sstr2_denovo_binder.py 2>&1

echo ""
echo "========================================="
echo "  All Arms Complete!"
echo "========================================="
