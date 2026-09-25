#!/bin/bash
#SBATCH --job-name=p1_fig5_mod_gdf15
#SBATCH --account=def-srussell
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=/scratch/srussell/srussell/logs/p1_fig5_mod_gdf15_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/p1_fig5_mod_gdf15_%j.err

set -euo pipefail

MANIFEST_DIR="/project/def-srussell/nirosha_xeno/scorial/01_state/slurm_manifests"
mkdir -p "$MANIFEST_DIR"
MANIFEST="$MANIFEST_DIR/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.yaml"

cat > "$MANIFEST" <<EOF
job_id: ${SLURM_JOB_ID}
job_name: ${SLURM_JOB_NAME}
task: P1 — re-derive Fig 5C/5D module-GDF15 rhos
script: /scratch/srussell/srussell/scripts/fig5_module_gdf15_rho.R
submitted: $(date -Is)
node: ${SLURMD_NODENAME:-unknown}
inputs:
  - /project/def-srussell/nirosha_xeno/wgcna/wgcna_BRAIN.rds
  - /project/def-srussell/nirosha_xeno/wgcna/module_assignments_BRAIN.tsv
  - /project/def-srussell/nirosha_xeno/deg/vsd_blind.rds
  - /project/def-srussell/nirosha_xeno/deg/sample_metadata.tsv
expected_outputs:
  - /project/def-srussell/nirosha_xeno/wgcna/brain_module_hallmark_enrichment.tsv
  - /project/def-srussell/nirosha_xeno/wgcna/brain_module_eigengene_trait_lm.tsv
  - /project/def-srussell/nirosha_xeno/wgcna/brain_module_gdf15L_correlations.tsv
  - /project/def-srussell/nirosha_xeno/figures/fig5/C_panel.pdf
  - /project/def-srussell/nirosha_xeno/figures/fig5/D_panel.pdf
  - /project/def-srussell/nirosha_xeno/figures/fig5/panel_metadata_C.yaml
  - /project/def-srussell/nirosha_xeno/figures/fig5/panel_metadata_D.yaml
EOF

module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21
export R_LIBS_USER="${HOME}/Rlibs/4.5:${R_LIBS_USER:-}"

cd /project/def-srussell/nirosha_xeno

Rscript /scratch/srussell/srussell/scripts/fig5_module_gdf15_rho.R

EXIT_CODE=$?
echo "exit_code: ${EXIT_CODE}" >> "$MANIFEST"
echo "finished: $(date -Is)" >> "$MANIFEST"
exit ${EXIT_CODE}
