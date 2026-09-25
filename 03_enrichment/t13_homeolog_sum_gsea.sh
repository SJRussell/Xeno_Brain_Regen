#!/bin/bash
#SBATCH --job-name=t13_homeolog_sum_gsea
#SBATCH --account=def-srussell
#SBATCH --time=01:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=2
#SBATCH --output=/scratch/srussell/srussell/logs/t13_homeolog_gsea_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/t13_homeolog_gsea_%j.err

set -euo pipefail
module purge
module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21
export R_LIBS_USER=/home/srussell/Rlibs/4.5

cd /project/def-srussell/nirosha_xeno
Rscript /scratch/srussell/srussell/scripts/t13_homeolog_sum_gsea.R

PROJ="/project/def-srussell/nirosha_xeno"
MANIFEST="$PROJ/scorial/01_state/slurm_manifests/job_${SLURM_JOB_ID}.yaml"
cat > "$MANIFEST" <<EOF
---
job_id: ${SLURM_JOB_ID}
job_name: t13_homeolog_sum_gsea
task_id: AUDIT-2026-05-05
stage: strand_homeolog_DE_audit
purpose: "GSEA on homeolog-summed primary DE results: MSigDB Hallmark + GO BP across 5 BRAIN + 9 BLAST contrasts."
end_time: $(date -Iseconds)
node: $(hostname)
exit_code: 0
outputs:
  - pathway/homeolog_sum_strand_corrected_enrichment.tsv
  - pathway/homeolog_sum_strand_corrected_results.rds
EOF
echo "Manifest: $MANIFEST"
