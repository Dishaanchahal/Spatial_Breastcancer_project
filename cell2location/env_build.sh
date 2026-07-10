#!/bin/bash
#SBATCH --job-name=c2l_env
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --output=/scratch/chahal.d/Spatial_GenAI/cell2location/logs/env_%j.out
#SBATCH --error=/scratch/chahal.d/Spatial_GenAI/cell2location/logs/env_%j.err

set -euo pipefail
module load miniconda3/25.9.1
source "$(conda info --base)/etc/profile.d/conda.sh"

ENV=/scratch/chahal.d/Spatial_GenAI/cell2location/c2l_env
cd /scratch/chahal.d/Spatial_GenAI/cell2location

if [ ! -d "$ENV" ]; then
  conda create -y -p "$ENV" -c conda-forge --override-channels python=3.10
fi
conda activate "$ENV"

pip install --no-cache-dir "cell2location[tutorials]" scanpy
python - <<'PY'
import cell2location, scanpy, torch, scvi
print("cell2location", cell2location.__version__)
print("scanpy", scanpy.__version__)
print("scvi-tools", scvi.__version__)
print("torch", torch.__version__)
PY
echo ENV_BUILD_DONE
