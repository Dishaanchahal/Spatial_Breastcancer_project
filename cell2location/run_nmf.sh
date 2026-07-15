#!/bin/bash
#SBATCH --job-name=c2l_nmf
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/chahal.d/Spatial_GenAI/cell2location/logs/nmf_%j.out
#SBATCH --error=/scratch/chahal.d/Spatial_GenAI/cell2location/logs/nmf_%j.err
set -euo pipefail
module load miniconda3/25.9.1
cd /scratch/chahal.d/Spatial_GenAI/cell2location
python3 17_nmf_microenvironments.py
echo NMF_PIPELINE_COMPLETE
