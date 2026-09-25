#!/bin/bash
#SBATCH --job-name=p2_fig1D_brain_pca
#SBATCH --account=def-srussell
#SBATCH --time=00:10:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/srussell/srussell/logs/p2_fig1D_brain_pca_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/p2_fig1D_brain_pca_%j.err

set -euo pipefail
MANIFEST="/project/def-srussell/nirosha_xeno/scorial/01_state/slurm_manifests/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.yaml"
cat > "$MANIFEST" <<EOF
job_id: ${SLURM_JOB_ID}
job_name: ${SLURM_JOB_NAME}
task: Fig 1.D — BRAIN PCA (homeolog-summed voom-normalized) PC1 vs PC2
script: /scratch/srussell/srussell/scripts/fig1D_brain_pca.R
submitted: $(date -Is)
expected_outputs:
  - /project/def-srussell/nirosha_xeno/figures/fig1/D_panel.pdf
  - /project/def-srussell/nirosha_xeno/figures/fig1/D_panel.png
  - /project/def-srussell/nirosha_xeno/figures/fig1/panel_metadata_D.yaml
  - /project/def-srussell/nirosha_xeno/deg/pca_BRAIN_homeolog_sum_strand_corrected.tsv
EOF

module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21
export R_LIBS_USER="${HOME}/Rlibs/4.5:${R_LIBS_USER:-}"
cd /project/def-srussell/nirosha_xeno
Rscript /scratch/srussell/srussell/scripts/fig1D_brain_pca.R

EXIT=$?
echo "exit_code: ${EXIT}"   >> "$MANIFEST"
echo "finished: $(date -Is)" >> "$MANIFEST"
exit ${EXIT}
