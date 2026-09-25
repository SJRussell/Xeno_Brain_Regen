#!/bin/bash
#SBATCH --job-name=t18_v10_brain_panels
#SBATCH --account=def-srussell
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=/scratch/srussell/srussell/logs/t18_v10_brain_panels_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/t18_v10_brain_panels_%j.err

set -euo pipefail
module purge
module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21
export R_LIBS_USER=/home/srussell/Rlibs/4.5

echo "[$(date -Iseconds)] t18 start (job ${SLURM_JOB_ID})"
Rscript /scratch/srussell/srussell/scripts/t18_v10_brain_panels.R
EXIT=$?
echo "[$(date -Iseconds)] t18 done (exit ${EXIT})"
exit ${EXIT}
