#!/bin/bash
#SBATCH --job-name=hA_walker_seurat
#SBATCH --account=def-srussell
#SBATCH --time=06:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --output=/scratch/srussell/srussell/logs/hA_walker_seurat_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/hA_walker_seurat_%j.err

# H-A Walker post-alevin: Seurat 5.x integrate + cluster + annotate + pseudobulk + UCell.
# Submit with --dependency=afterok:<alevin_array_job_id> to auto-launch.

set -euo pipefail
module purge
module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21
export R_LIBS_USER=/home/srussell/Rlibs/4.5

cd /project/def-srussell/nirosha_xeno
Rscript /scratch/srussell/srussell/scripts/hA_walker_seurat.R

PROJ="/project/def-srussell/nirosha_xeno"
MANIFEST="$PROJ/scorial/01_state/slurm_manifests/job_${SLURM_JOB_ID}.yaml"
cat > "$MANIFEST" <<EOF
---
job_id: ${SLURM_JOB_ID}
job_name: hA_walker_seurat
task_id: T-09a-seurat
stage: comparator_walker_seurat
purpose: "DEC-017 steps 4-7: Seurat 5.x integrate 4 axolotl brain samples; cluster + Walker-marker annotation; pseudobulk; UCell-score 4 Xenopus brain CT signatures + GDF15 against axolotl populations; treatment-population enrichment for etv1+ cluster."
end_time: $(date -Iseconds)
node: $(hostname)
exit_code: 0
outputs:
  - comparators/walker/walker_integrated.rds
  - comparators/walker/cluster_annotations.tsv
  - comparators/walker/cluster_markers.tsv
  - comparators/walker/pseudobulk_counts.tsv
  - comparators/walker/ucell_signatures_per_cell.tsv
  - comparators/walker/ucell_signatures_per_cluster.tsv
  - comparators/walker/treatment_pop_enrichment.tsv
  - figures/stage_H/walker/walker_umap_clusters.pdf
  - figures/stage_H/walker/walker_umap_treatment.pdf
  - figures/stage_H/walker/walker_signature_heatmap.pdf
EOF
echo "Manifest: $MANIFEST"
