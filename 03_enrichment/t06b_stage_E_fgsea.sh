#!/bin/bash
#SBATCH --job-name=t06b_fgsea
#SBATCH --account=def-srussell
#SBATCH --time=01:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --output=/scratch/srussell/srussell/logs/t06b_fgsea_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/t06b_fgsea_%j.err

set -euo pipefail
module purge
module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21
export R_LIBS_USER=/home/srussell/Rlibs/4.5

cd /project/def-srussell/nirosha_xeno
echo "[$(date -Iseconds)] t06b start (job ${SLURM_JOB_ID})"
Rscript /scratch/srussell/srussell/scripts/t06b_stage_E_fgsea.R
EXIT=$?
echo "[$(date -Iseconds)] t06b done (exit ${EXIT})"

MAN=/project/def-srussell/nirosha_xeno/scorial/01_state/slurm_manifests/job_${SLURM_JOB_ID}.yaml
cat > "$MAN" <<EOF
---
job_id: ${SLURM_JOB_ID}
task_id: T-06b
stage: fgsea_enrichment_deseq2_primary
exit_code: ${EXIT}
end_time: $(date -Iseconds)
tool: fgsea (r-bundle-bioconductor/3.21)
inputs:
  - deg/deseq2_BRAIN_primary.rds
  - deg/deseq2_BLAST_primary.rds
  - deg/voom_BRAIN_primary.rds  (confirmatory)
  - deg/voom_BLAST_primary.rds  (confirmatory)
outputs:
  - pathway/*fgsea*.xlsx
  - pathway/all_fgsea_summary.tsv
EOF
exit ${EXIT}
