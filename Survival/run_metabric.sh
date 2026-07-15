#!/bin/bash
#SBATCH --job-name=metabric_surv
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:40:00
#SBATCH --output=/scratch/chahal.d/Spatial_GenAI/Survival/logs/metabric_%j.out
#SBATCH --error=/scratch/chahal.d/Spatial_GenAI/Survival/logs/metabric_%j.err
set -euo pipefail
ENV=/scratch/chahal.d/Spatial_GenAI/CellComm/ccc_env
export LD_LIBRARY_PATH="$ENV/lib:${LD_LIBRARY_PATH:-}"
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
"$ENV/bin/python" -c "import bz2; print('bz2 import OK')"
cd /scratch/chahal.d/Spatial_GenAI/Survival
"$ENV/bin/python" 16_survival_metabric.py
echo METABRIC_PIPELINE_COMPLETE
