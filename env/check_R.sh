#!/bin/bash
#SBATCH --job-name=check_R
#SBATCH --account=def-srussell
#SBATCH --time=00:10:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/srussell/srussell/logs/check_R_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/check_R_%j.err
set -euo pipefail
module purge
module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21
Rscript /scratch/srussell/srussell/scripts/check_R_packages.R
