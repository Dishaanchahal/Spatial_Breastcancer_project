#!/bin/bash
#SBATCH --job-name=spacexr_install3
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --output=/scratch/chahal.d/Spatial_GenAI/Validation_RCTD/logs/install3_%j.out
#SBATCH --error=/scratch/chahal.d/Spatial_GenAI/Validation_RCTD/logs/install3_%j.err
set -euo pipefail
module load R/4.4.1
module load OpenBLAS/0.3.29
export LD_LIBRARY_PATH=$HOME/rdeps/lib:$LD_LIBRARY_PATH
export R_MAKEVARS_USER=$HOME/rdeps/Makevars
export MAKEFLAGS="-j8"
# purge stale binaries built with the wrong toolchain, then rebuild clean
Rscript -e 'for (p in c("spacexr","fields","spam","spam64")) try(remove.packages(p), silent=TRUE)'
Rscript -e 'options(Ncpus=8); install.packages(c("spam","fields"), repos="https://cloud.r-project.org")'
Rscript -e 'suppressPackageStartupMessages({library(spam); library(fields)}); cat("FIELDS_OK\n")'
Rscript -e 'options(Ncpus=8); remotes::install_github("dmcable/spacexr", upgrade="never", build_vignettes=FALSE)'
Rscript -e 'suppressPackageStartupMessages(library(spacexr)); cat("SPACEXR_OK\n")'
echo INSTALL3_DONE
