#!/bin/bash
#SBATCH --job-name=t05c_deseq2_wald
#SBATCH --account=def-srussell
#SBATCH --time=00:45:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=/scratch/srussell/srussell/logs/t05c_deseq2_wald_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/t05c_deseq2_wald_%j.err

set -euo pipefail
module purge
module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21
export R_LIBS_USER=/home/srussell/Rlibs/4.5

MAN=/project/def-srussell/nirosha_xeno/scorial/01_state/slurm_manifests
mkdir -p "$MAN"

echo "[$(date -Iseconds)] t05c start (job ${SLURM_JOB_ID})"
Rscript /scratch/srussell/srussell/scripts/t05c_deseq2_wald.R
EXIT=$?
echo "[$(date -Iseconds)] t05c done (exit ${EXIT})"

cat > "${MAN}/job_${SLURM_JOB_ID}.yaml" <<EOF
---
job_id: ${SLURM_JOB_ID}
task_id: T-05c
stage: deseq2_wald_primary_DE
exit_code: ${EXIT}
end_time: $(date -Iseconds)
tool: DESeq2 (r-bundle-bioconductor/3.21)
parameters:
  design: '~ 0 + grp  (grp = treatment_time per tissue)'
  test: Wald
  prefilter: 'CPM > 1 in >= 3 samples (matches limma pipeline filter)'
  independent_filtering: TRUE
  alpha: 0.05
outputs:
  - deg/deseq2_BRAIN_primary.{rds,xlsx,_summary.tsv}
  - deg/deseq2_BLAST_primary.{rds,xlsx,_summary.tsv}
  - deg/deseq2_BRAIN_sensitivity.{rds,xlsx,_summary.tsv}
  - deg/all_deseq2_wald_summary.tsv
EOF
exit ${EXIT}
