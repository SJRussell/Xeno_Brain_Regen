#!/bin/bash
#SBATCH --job-name=t10_homeolog_sum
#SBATCH --account=def-srussell
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=/scratch/srussell/srussell/logs/t10_homeolog_sum_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/t10_homeolog_sum_%j.err

set -euo pipefail
module purge
module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21
export R_LIBS_USER=/home/srussell/Rlibs/4.5

cd /project/def-srussell/nirosha_xeno
Rscript /scratch/srussell/srussell/scripts/t10_homeolog_sum_de.R

PROJ="/project/def-srussell/nirosha_xeno"
MANIFEST="$PROJ/scorial/01_state/slurm_manifests/job_${SLURM_JOB_ID}.yaml"
cat > "$MANIFEST" <<EOF
---
job_id: ${SLURM_JOB_ID}
job_name: t10_homeolog_sum
task_id: T-10
stage: stage_D_sensitivity_homeolog
purpose: "Homeolog-summed sensitivity DE per Stewart's 2026-05-01 spec. SECONDARY analysis (does NOT overwrite primary). Sums .L/.S featureCounts before limma-voom; numbered paralogs (X.n.L/S) summed separately under X.n. Outputs voom_{BRAIN,BLAST}_homeolog_sum.{xlsx,rds} and v1 manuscript concordance table."
end_time: $(date -Iseconds)
node: $(hostname)
exit_code: 0
outputs:
  - deg/voom_BRAIN_homeolog_sum.xlsx
  - deg/voom_BLAST_homeolog_sum.xlsx
  - deg/all_de_summary_homeolog_sum.tsv
  - deg/v1_manuscript_concordance_homeolog_sum.tsv
EOF
echo "Manifest: $MANIFEST"
