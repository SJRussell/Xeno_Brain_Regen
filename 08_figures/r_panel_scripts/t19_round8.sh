#!/bin/bash
#SBATCH --job-name=t19_round8
#SBATCH --account=def-srussell
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=/scratch/srussell/srussell/logs/t19_round8_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/t19_round8_%j.err

set -euo pipefail
module purge
module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21
export R_LIBS_USER=/home/srussell/Rlibs/4.5

echo "[$(date -Iseconds)] t19 start (job ${SLURM_JOB_ID})"
Rscript /scratch/srussell/srussell/scripts/t19_round8.R
EXIT=$?
echo "[$(date -Iseconds)] t19 done (exit ${EXIT})"
exit ${EXIT}
