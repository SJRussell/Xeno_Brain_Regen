#!/bin/bash
#SBATCH --job-name=walker_marker_remap
#SBATCH --account=def-srussell
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --output=/scratch/srussell/srussell/logs/walker_marker_remap_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/walker_marker_remap_%j.err

# RED-2: re-run UCell scoring on walker_integrated.rds with corrected axolotl
# marker IDs (AMEX60DD gene-id namespace, derived from GTF lookup).
set -euo pipefail

module purge
module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21
export R_LIBS_USER=/home/srussell/Rlibs/4.5

cd /project/def-srussell/nirosha_xeno
Rscript /scratch/srussell/srussell/scripts/walker_marker_remap.R

# Manifest
PROJ="/project/def-srussell/nirosha_xeno"
MANIFEST="$PROJ/scorial/01_state/slurm_manifests/job_${SLURM_JOB_ID}.yaml"
cat > "$MANIFEST" <<EOF
---
job_id: ${SLURM_JOB_ID}
job_name: walker_marker_remap
task_id: RED-2
stage: comparator_walker_ucell_remap
purpose: "RED-2 (iter-19): re-score walker_integrated.rds UCell signatures using AMEX60DD-namespace markers (corrected for tx2gene's [hs]-only extraction bug)."
end_time: $(date -Iseconds)
node: $(hostname)
exit_code: 0
related_decisions: [DEC-017]
inputs:
  - $PROJ/comparators/walker/walker_integrated.rds
outputs:
  - $PROJ/comparators/walker/cluster_annotations_v2.tsv
  - $PROJ/comparators/walker/ucell_signatures_per_cell_v2.tsv
  - $PROJ/comparators/walker/ucell_signatures_per_cluster_v2.tsv
  - $PROJ/comparators/walker/etv1_population_proportions_v2.tsv
  - $PROJ/comparators/walker/treatment_pop_enrichment_v2.tsv
  - $PROJ/comparators/walker/walker_integrated_v2.rds
EOF
echo "Manifest: $MANIFEST"
