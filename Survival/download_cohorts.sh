#!/bin/bash
#SBATCH --job-name=surv_dl
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:40:00
#SBATCH --output=/scratch/chahal.d/Spatial_GenAI/Survival/logs/dl_%j.out
#SBATCH --error=/scratch/chahal.d/Spatial_GenAI/Survival/logs/dl_%j.err
set -euo pipefail
cd /scratch/chahal.d/Spatial_GenAI/Survival
for d in brca_metabric brca_tcga_pan_can_atlas_2018; do
  curl -sL -o $d.tar.gz "https://cbioportal-datahub.s3.amazonaws.com/$d.tar.gz"
  tar xzf $d.tar.gz && rm $d.tar.gz
  echo "$d files:"; ls $d | grep -iE "mrna|clinical" | head
done
echo SURV_DOWNLOAD_DONE
