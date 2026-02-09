#!/bin/bash
# Test MolMIM hosted API with molmim.key (nvapi- format)
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate bio-tools 2>/dev/null

# molmim.key 우선, 없으면 ngc.key
if [ -f /mnt/g/repos/bio/molmim.key ]; then
  API_KEY=$(cat /mnt/g/repos/bio/molmim.key)
  echo "Using: molmim.key"
elif [ -f /mnt/g/repos/bio/ngc.key ]; then
  API_KEY=$(cat /mnt/g/repos/bio/ngc.key)
  echo "Using: ngc.key"
else
  echo "ERROR: No key file found"
  exit 1
fi
echo "Key prefix: ${API_KEY:0:10}..."

echo ""
echo "=== Test: /generate endpoint ==="
curl -s -w "\nHTTP_STATUS: %{http_code}\n" \
  -X POST \
  "https://health.api.nvidia.com/v1/biology/nvidia/molmim/generate" \
  -H "accept: application/json" \
  -H "Content-Type: application/json" \
  -H "Authorization: Bearer ${API_KEY}" \
  -d '{"smi": "CCO", "algorithm": "none", "num_molecules": 3, "particles": 8, "scaled_radius": 1.0}' 2>&1
