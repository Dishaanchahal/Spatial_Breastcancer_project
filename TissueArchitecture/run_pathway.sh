#!/bin/bash
#SBATCH --job-name=pathway
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/chahal.d/Spatial_GenAI/TissueArchitecture/logs/pathway_%j.out
#SBATCH --error=/scratch/chahal.d/Spatial_GenAI/TissueArchitecture/logs/pathway_%j.err
set -euo pipefail
cd /scratch/chahal.d/Spatial_GenAI/TissueArchitecture
/scratch/chahal.d/Spatial_GenAI/CellComm/ccc_env/bin/python 12_pathway_activity.py
echo PATHWAY_PIPELINE_COMPLETE
