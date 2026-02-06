#!/usr/bin/env bash
# Install Miniconda (Linux) inside WSL.
# Run from WSL: bash scripts/install_miniconda_wsl.sh
# Then run: source ~/.bashrc  (or open a new terminal)

set -e
INSTALLER="Miniconda3-latest-Linux-x86_64.sh"
URL="https://repo.anaconda.com/miniconda/${INSTALLER}"

echo "Downloading Miniconda for Linux..."
wget -q "$URL" -O "$INSTALLER"

echo "Installing Miniconda (use default prefix, say yes to init)..."
bash "$INSTALLER" -b -p "$HOME/miniconda3"
rm -f "$INSTALLER"

# Init so conda is on PATH in new shells
"$HOME/miniconda3/bin/conda" init bash

echo "Done. Run: source ~/.bashrc   then: conda --version"
