#!/bin/bash
#SBATCH --job-name=he_within
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/chahal.d/Spatial_GenAI/HE_FoundationModel/logs/within_%j.out
#SBATCH --error=/scratch/chahal.d/Spatial_GenAI/HE_FoundationModel/logs/within_%j.err
set -euo pipefail
module load miniconda3/25.9.1
cd /scratch/chahal.d/Spatial_GenAI/HE_FoundationModel
python3 he_predict_within.py
