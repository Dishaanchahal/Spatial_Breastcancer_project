#!/bin/bash
#SBATCH --job-name=tissue_arch
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/chahal.d/Spatial_GenAI/TissueArchitecture/logs/tissue_%j.out
#SBATCH --error=/scratch/chahal.d/Spatial_GenAI/TissueArchitecture/logs/tissue_%j.err
set -euo pipefail
cd /scratch/chahal.d/Spatial_GenAI/TissueArchitecture
PY=/scratch/chahal.d/Spatial_GenAI/CellComm/ccc_env/bin/python
echo "[$(date)] tissue architecture"; "$PY" 11_tissue_architecture.py
echo "[$(date)] pathway activity";   "$PY" 12_pathway_activity.py
echo "[$(date)] TISSUE_PIPELINE_COMPLETE"
