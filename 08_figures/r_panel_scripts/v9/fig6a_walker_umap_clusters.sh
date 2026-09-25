#!/bin/bash
#SBATCH --job-name=v9_fig6a
#SBATCH --account=def-srussell
#SBATCH --time=00:30:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=2
#SBATCH --output=/scratch/srussell/srussell/logs/v9_fig6a_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/v9_fig6a_%j.err

set -euo pipefail
MANIFEST="/project/def-srussell/nirosha_xeno/scorial/01_state/slurm_manifests/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.yaml"
cat > "$MANIFEST" <<EOF
job_id: ${SLURM_JOB_ID}
job_name: ${SLURM_JOB_NAME}
task: Fig 6.a (v9) - Walker axolotl brain scRNA UMAP, 23 clusters (d3 cat20 palette extended)
script: /scratch/srussell/srussell/scripts/v9/fig6a_walker_umap_clusters.R
submitted: $(date -Is)
expected_outputs:
  - /project/def-srussell/nirosha_xeno/figures_v9/fig6/a_panel.pdf
  - /project/def-srussell/nirosha_xeno/figures_v9/fig6/a_panel.png
  - /project/def-srussell/nirosha_xeno/figures_v9/fig6/panel_meta_a.yaml
EOF

module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21
export R_LIBS_USER="${HOME}/Rlibs/4.5:${R_LIBS_USER:-}"
cd /project/def-srussell/nirosha_xeno
Rscript /scratch/srussell/srussell/scripts/v9/fig6a_walker_umap_clusters.R

EXIT=$?
echo "exit_code: ${EXIT}"   >> "$MANIFEST"
echo "finished: $(date -Is)" >> "$MANIFEST"
exit ${EXIT}
