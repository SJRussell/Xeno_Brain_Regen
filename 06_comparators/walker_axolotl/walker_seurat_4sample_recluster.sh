#!/bin/bash
#SBATCH --job-name=walker_4sample
#SBATCH --account=def-srussell
#SBATCH --time=06:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --output=/scratch/srussell/srussell/logs/walker_4sample_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/walker_4sample_%j.err

# Re-run Seurat integration with all 4 samples (post-RED-1 muscimol recovery).
# Submitted with --dependency=afterok:<alevin-musc job>.
# This re-uses the existing hA_walker_seurat.R but operates on 4-sample alevin
# output (3-sample existing + new muscimol).
set -euo pipefail
module purge
module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21
export R_LIBS_USER=/home/srussell/Rlibs/4.5

# Backup the prior 3-sample integration object before overwriting
PROJ="/project/def-srussell/nirosha_xeno"
PROJ_OUT="$PROJ/comparators/walker"
if [ -f "$PROJ_OUT/walker_integrated.rds" ]; then
    BAK="$PROJ_OUT/walker_integrated_3sample_$(date +%Y%m%dT%H%M%S).rds"
    cp "$PROJ_OUT/walker_integrated.rds" "$BAK"
    echo "[$(date -Iseconds)] preserved prior 3-sample integration as $BAK"
fi

cd /project/def-srussell/nirosha_xeno
Rscript /scratch/srussell/srussell/scripts/hA_walker_seurat.R

MANIFEST="$PROJ/scorial/01_state/slurm_manifests/job_${SLURM_JOB_ID}.yaml"
cat > "$MANIFEST" <<EOF
---
job_id: ${SLURM_JOB_ID}
job_name: walker_4sample
task_id: RED-1-seurat-4sample
stage: comparator_walker_seurat_4sample_recluster
purpose: "RED-1 follow-up: re-run Seurat 5 SCT/CCA integration with all 4 Walker samples (uninjured + tail + limb + muscimol) after RED-1 recovered SRR31800874. Same pipeline params as hA_walker_seurat.R (vst.flavor=v1, sequential plan, FindClusters resolution 0.6). Re-derives walker_integrated.rds with muscimol cells included."
end_time: $(date -Iseconds)
node: $(hostname)
exit_code: 0
related_decisions: [DEC-017, DEC-019]
outputs:
  - comparators/walker/walker_integrated.rds
  - comparators/walker/cluster_annotations.tsv
  - comparators/walker/cluster_markers.tsv
  - comparators/walker/pseudobulk_counts.tsv
  - comparators/walker/ucell_signatures_per_cell.tsv
  - comparators/walker/ucell_signatures_per_cluster.tsv
  - comparators/walker/treatment_pop_enrichment.tsv
EOF
echo "Manifest: $MANIFEST"
