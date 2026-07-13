#!/bin/bash
#SBATCH --job-name=ccc_env
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --output=/scratch/chahal.d/Spatial_GenAI/CellComm/logs/env_%j.out
#SBATCH --error=/scratch/chahal.d/Spatial_GenAI/CellComm/logs/env_%j.err
set -euo pipefail
module load miniconda3/25.9.1
source "$(conda info --base)/etc/profile.d/conda.sh"
ENV=/scratch/chahal.d/Spatial_GenAI/CellComm/ccc_env
[ -d "$ENV" ] || conda create -y -p "$ENV" -c conda-forge --override-channels python=3.10
"$ENV/bin/pip" install --no-cache-dir liana scanpy
"$ENV/bin/python" - <<'PY'
import liana, scanpy, pandas, numpy
print("liana", liana.__version__, "| scanpy", scanpy.__version__, "| pandas", pandas.__version__)
PY
echo CCC_ENV_DONE
