#!/bin/bash
# ============================================================================
# run_pipeline.sh  —  SLURM driver for the spatial transcriptomics pipeline
#
# Submit from the project root on Explorer:
#     cd /scratch/chahal.d/Spatial_GenAI
#     sbatch run_pipeline.sh
#
# Or run a single step interactively:
#     module load R/4.4.1 && Rscript 03_label_transfer.R
# ============================================================================
#SBATCH --job-name=spatial_brca
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/spatial_%j.out
#SBATCH --error=logs/spatial_%j.err

set -euo pipefail

module load R/4.4.1

cd /scratch/chahal.d/Spatial_GenAI
mkdir -p logs

echo "[$(date)] Step 1: process spatial samples"
Rscript 01_process_spatial.R

echo "[$(date)] Step 2: build scRNA-seq reference"
Rscript 02_process_reference.R

echo "[$(date)] Step 3: label transfer"
Rscript 03_label_transfer.R

echo "[$(date)] Step 4: label-transfer plots"
Rscript 04_label_transfer_plots.R

echo "[$(date)] Step 5: pathologist annotations"
Rscript 05_pathologist_annotations.R

echo "[$(date)] Step 6: validate annotations"
Rscript 06_validate_annotations.R

echo "[$(date)] Step 7: TME comparison"
Rscript 07_tme_comparison.R

echo "[$(date)] Pipeline complete."
