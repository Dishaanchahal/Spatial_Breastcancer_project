#!/bin/bash
#SBATCH --job-name=he_fm
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/chahal.d/Spatial_GenAI/HE_FoundationModel/logs/he_%j.out
#SBATCH --error=/scratch/chahal.d/Spatial_GenAI/HE_FoundationModel/logs/he_%j.err
set -euo pipefail
module load miniconda3/25.9.1
export HF_HOME=/scratch/chahal.d/hf_cache
cd /scratch/chahal.d/Spatial_GenAI/HE_FoundationModel
echo "[$(date)] embedding H&E patches with Phikon"; python3 he_embed.py
echo "[$(date)] predicting composition + TLS from morphology"; python3 he_predict.py
echo "[$(date)] HE_PIPELINE_COMPLETE"
