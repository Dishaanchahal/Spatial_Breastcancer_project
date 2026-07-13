#!/bin/bash
#SBATCH --job-name=survival
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:40:00
#SBATCH --output=/scratch/chahal.d/Spatial_GenAI/Survival/logs/surv_%j.out
#SBATCH --error=/scratch/chahal.d/Spatial_GenAI/Survival/logs/surv_%j.err
set -euo pipefail
cd /scratch/chahal.d/Spatial_GenAI/Survival
/scratch/chahal.d/Spatial_GenAI/CellComm/ccc_env/bin/python 13_survival_signatures.py
echo SURVIVAL_PIPELINE_COMPLETE
