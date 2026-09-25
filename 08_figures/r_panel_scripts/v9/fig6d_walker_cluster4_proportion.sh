#!/bin/bash
#SBATCH --job-name=v9_fig6d
#SBATCH --account=def-srussell
#SBATCH --time=00:30:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=2
#SBATCH --output=/scratch/srussell/srussell/logs/v9_fig6d_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/v9_fig6d_%j.err

set -euo pipefail
MANIFEST="/project/def-srussell/nirosha_xeno/scorial/01_state/slurm_manifests/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.yaml"
cat > "$MANIFEST" <<EOF
job_id: ${SLURM_JOB_ID}
job_name: ${SLURM_JOB_NAME}
task: Fig 6.d (v9) - Cluster 4 proportion across 4 Walker treatments + Fisher OR vs uninjured
script: /scratch/srussell/srussell/scripts/v9/fig6d_walker_cluster4_proportion.R
submitted: $(date -Is)
expected_outputs:
  - /project/def-srussell/nirosha_xeno/figures_v9/fig6/d_panel.pdf
  - /project/def-srussell/nirosha_xeno/figures_v9/fig6/d_panel.png
  - /project/def-srussell/nirosha_xeno/figures_v9/fig6/panel_meta_d.yaml
  - /project/def-srussell/nirosha_xeno/figures_v9/fig6/fig6d_cluster4_proportion_fisher.tsv
EOF

module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21
export R_LIBS_USER="${HOME}/Rlibs/4.5:${R_LIBS_USER:-}"
cd /project/def-srussell/nirosha_xeno
Rscript /scratch/srussell/srussell/scripts/v9/fig6d_walker_cluster4_proportion.R

EXIT=$?
echo "exit_code: ${EXIT}"   >> "$MANIFEST"
echo "finished: $(date -Is)" >> "$MANIFEST"
exit ${EXIT}
