#!/bin/bash
#SBATCH --job-name=c2l_fix
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:40:00
#SBATCH --output=/scratch/chahal.d/Spatial_GenAI/cell2location/logs/fix_%j.out
#SBATCH --error=/scratch/chahal.d/Spatial_GenAI/cell2location/logs/fix_%j.err
set -euo pipefail
module load miniconda3/25.9.1
echo "before:"; python3 -c "import torch;print(torch.__version__, torch.version.cuda)" 2>&1 | tail -1
python3 -m pip install --user --no-cache-dir --force-reinstall torch --index-url https://download.pytorch.org/whl/cu124
echo "after:"; python3 -c "import torch;print('torch', torch.__version__, 'cuda', torch.version.cuda)"
echo FIX_DONE
