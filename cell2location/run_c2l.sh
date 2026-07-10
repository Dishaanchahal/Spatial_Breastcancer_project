#!/bin/bash
#SBATCH --job-name=c2l_run
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/chahal.d/Spatial_GenAI/cell2location/logs/run_%j.out
#SBATCH --error=/scratch/chahal.d/Spatial_GenAI/cell2location/logs/run_%j.err

set -euo pipefail
module load miniconda3/25.9.1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /scratch/chahal.d/Spatial_GenAI/cell2location/c2l_env

cd /scratch/chahal.d/Spatial_GenAI/cell2location

echo "[$(date)] Stage 0: prep AnnData"
python c2l_prep.py

echo "[$(date)] Stage 1+2: reference regression + spatial mapping (GPU)"
python c2l_run.py

echo "[$(date)] Stage 3: plots + composition"
python c2l_plots.py

echo "[$(date)] C2L_PIPELINE_COMPLETE"
