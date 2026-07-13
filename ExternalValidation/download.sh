#!/bin/bash
#SBATCH --job-name=xval_dl
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/chahal.d/Spatial_GenAI/ExternalValidation/logs/dl_%j.out
#SBATCH --error=/scratch/chahal.d/Spatial_GenAI/ExternalValidation/logs/dl_%j.err
set -euo pipefail
base=https://cf.10xgenomics.com/samples/spatial-exp/1.1.0
cd /scratch/chahal.d/Spatial_GenAI/ExternalValidation/samples
for s in V1_Breast_Cancer_Block_A_Section_1 V1_Breast_Cancer_Block_A_Section_2; do
  mkdir -p "$s"; cd "$s"
  curl -sL -o m.tar.gz  "$base/$s/${s}_filtered_feature_bc_matrix.tar.gz"
  curl -sL -o sp.tar.gz "$base/$s/${s}_spatial.tar.gz"
  tar xzf m.tar.gz && tar xzf sp.tar.gz && rm m.tar.gz sp.tar.gz
  echo "$s:"; ls filtered_feature_bc_matrix spatial
  cd ..
done
echo XVAL_DOWNLOAD_DONE
