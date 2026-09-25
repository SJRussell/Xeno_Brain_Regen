#!/bin/bash
#SBATCH --job-name=hD_liao_an
#SBATCH --account=def-srussell
#SBATCH --time=01:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=/scratch/srussell/srussell/logs/hD_liao_an_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/hD_liao_an_%j.err

set -euo pipefail
module purge
module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21
export R_LIBS_USER=/home/srussell/Rlibs/4.5

cd /project/def-srussell/nirosha_xeno
Rscript /scratch/srussell/srussell/scripts/hD_liao_analyze.R

MANIFEST="/project/def-srussell/nirosha_xeno/scorial/01_state/slurm_manifests/job_${SLURM_JOB_ID}.yaml"
cat > "$MANIFEST" <<EOF
---
job_id: ${SLURM_JOB_ID}
job_name: hD_liao_an
task_id: T-09d-analyze
stage: comparator_liao_analysis
purpose: "Liao 2022 brain atlas — per-celltype marker derivation (CPM-pseudobulk log-FC, top 50/celltype, ≥30 cells/celltype filter); UCell-on-bulk scoring vs Xenopus 12 brain samples; treatment×time interaction lm per celltype."
end_time: $(date -Iseconds)
node: $(hostname)
exit_code: 0
outputs:
  workspace: /scratch/srussell/srussell/nirosha_xeno_comparators/liao/work
  per_celltype_markers: /project/def-srussell/nirosha_xeno/comparators/liao/per_celltype_markers.tsv
  per_celltype_scores: /project/def-srussell/nirosha_xeno/comparators/liao/per_celltype_scores.tsv
  per_celltype_interaction: /project/def-srussell/nirosha_xeno/comparators/liao/per_celltype_interaction.tsv
EOF
echo "Manifest written: $MANIFEST"
