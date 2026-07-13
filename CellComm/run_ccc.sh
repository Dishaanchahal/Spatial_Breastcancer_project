#!/bin/bash
#SBATCH --job-name=ccc_run
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/chahal.d/Spatial_GenAI/CellComm/logs/run_%j.out
#SBATCH --error=/scratch/chahal.d/Spatial_GenAI/CellComm/logs/run_%j.err
set -euo pipefail
cd /scratch/chahal.d/Spatial_GenAI/CellComm
ENV=/scratch/chahal.d/Spatial_GenAI/CellComm/ccc_env
"$ENV/bin/python" 10_cell_communication.py
echo CCC_PIPELINE_COMPLETE
