#!/usr/bin/env bash
set -e
cd /mnt/g/repos/bio

echo "=== Checking gh CLI ==="
if ! command -v gh &> /dev/null; then
    echo "[INFO] gh not found. Installing..."
    type -p curl >/dev/null || sudo apt install curl -y
    curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg | sudo dd of=/usr/share/keyrings/githubcli-archive-keyring.gpg
    sudo chmod go+r /usr/share/keyrings/githubcli-archive-keyring.gpg
    echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" | sudo tee /etc/apt/sources.list.d/github-cli.list > /dev/null
    sudo apt update
    sudo apt install gh -y
fi

echo "=== gh version ==="
gh --version

echo "=== Checking auth status ==="
gh auth status || echo "[WARN] Not authenticated. Run: gh auth login"
