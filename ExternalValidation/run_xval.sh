#!/bin/bash
#SBATCH --job-name=xval_run
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/chahal.d/Spatial_GenAI/ExternalValidation/logs/xval_%j.out
#SBATCH --error=/scratch/chahal.d/Spatial_GenAI/ExternalValidation/logs/xval_%j.err
set -euo pipefail
module load miniconda3/25.9.1
module load OpenBLAS/0.3.29 2>/dev/null || true
export LD_LIBRARY_PATH=$HOME/rdeps/lib:$LD_LIBRARY_PATH
cd /scratch/chahal.d/Spatial_GenAI/ExternalValidation
echo "[$(date)] cell2location mapping (reuse reference signatures)"; python3 xval_c2l.py
echo "[$(date)] downstream (composition/TLS/PROGENy)"; /scratch/chahal.d/Spatial_GenAI/CellComm/ccc_env/bin/python xval_downstream.py
echo "[$(date)] XVAL_PIPELINE_COMPLETE"
