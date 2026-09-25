#!/bin/bash
#SBATCH --job-name=v9_fig6c
#SBATCH --account=def-srussell
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/srussell/srussell/logs/v9_fig6c_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/v9_fig6c_%j.err

set -euo pipefail
MANIFEST="/project/def-srussell/nirosha_xeno/scorial/01_state/slurm_manifests/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.yaml"
cat > "$MANIFEST" <<EOF
job_id: ${SLURM_JOB_ID}
job_name: ${SLURM_JOB_NAME}
task: Fig 6.c (v9) - Walker UCell etv1+/glut score by cluster, ranked bar
script: /scratch/srussell/srussell/scripts/v9/fig6c_walker_ucell_ranked_bar.R
submitted: $(date -Is)
expected_outputs:
  - /project/def-srussell/nirosha_xeno/figures_v9/fig6/c_panel.pdf
  - /project/def-srussell/nirosha_xeno/figures_v9/fig6/c_panel.png
  - /project/def-srussell/nirosha_xeno/figures_v9/fig6/panel_meta_c.yaml
EOF

module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21
export R_LIBS_USER="${HOME}/Rlibs/4.5:${R_LIBS_USER:-}"
cd /project/def-srussell/nirosha_xeno
Rscript /scratch/srussell/srussell/scripts/v9/fig6c_walker_ucell_ranked_bar.R

EXIT=$?
echo "exit_code: ${EXIT}"   >> "$MANIFEST"
echo "finished: $(date -Is)" >> "$MANIFEST"
exit ${EXIT}
