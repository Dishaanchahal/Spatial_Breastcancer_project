#!/bin/bash
#SBATCH --job-name=spacexr_install
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --output=/scratch/chahal.d/Spatial_GenAI/Validation_RCTD/logs/install_%j.out
#SBATCH --error=/scratch/chahal.d/Spatial_GenAI/Validation_RCTD/logs/install_%j.err
set -euo pipefail
module load R/4.4.1
Rscript -e 'if(!requireNamespace("remotes",quietly=TRUE)) install.packages("remotes",repos="https://cloud.r-project.org"); options(Ncpus=8); remotes::install_github("dmcable/spacexr", upgrade="never", build_vignettes=FALSE)'
Rscript -e 'suppressPackageStartupMessages(library(spacexr)); cat("SPACEXR_OK\n")'
echo INSTALL_DONE
