#!/bin/bash
#SBATCH --job-name=rctd_run
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/chahal.d/Spatial_GenAI/Validation_RCTD/logs/rctd_%j.out
#SBATCH --error=/scratch/chahal.d/Spatial_GenAI/Validation_RCTD/logs/rctd_%j.err
set -euo pipefail
module load R/4.4.1
module load OpenBLAS/0.3.29
export LD_LIBRARY_PATH=$HOME/rdeps/lib:$LD_LIBRARY_PATH
cd /scratch/chahal.d/Spatial_GenAI
echo "[$(date)] RCTD run"
Rscript 08_rctd_validation.R
echo "[$(date)] compare"
Rscript 09_compare_c2l_rctd.R
echo "[$(date)] RCTD_PIPELINE_COMPLETE"
