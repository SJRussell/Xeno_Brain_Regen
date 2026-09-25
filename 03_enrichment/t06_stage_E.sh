#!/bin/bash
#SBATCH --job-name=t06_stage_E
#SBATCH --account=def-srussell
#SBATCH --time=01:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=/scratch/srussell/srussell/logs/t06_stage_E_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/t06_stage_E_%j.err

set -euo pipefail
module purge
module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21

export R_LIBS_USER=/home/srussell/Rlibs/4.5

cd /project/def-srussell/nirosha_xeno
Rscript /scratch/srussell/srussell/scripts/t06_stage_E_enrichment.R
