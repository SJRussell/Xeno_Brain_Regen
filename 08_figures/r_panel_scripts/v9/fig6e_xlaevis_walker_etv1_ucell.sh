#!/bin/bash
#SBATCH --job-name=v9_fig6e
#SBATCH --account=def-srussell
#SBATCH --time=00:15:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/srussell/srussell/logs/v9_fig6e_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/v9_fig6e_%j.err

set -euo pipefail
MANIFEST="/project/def-srussell/nirosha_xeno/scorial/01_state/slurm_manifests/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.yaml"
cat > "$MANIFEST" <<EOF
job_id: ${SLURM_JOB_ID}
job_name: ${SLURM_JOB_NAME}
task: Fig 6.e (v9) - X. laevis BRAIN UCell vs Walker etv1+ signature, per-homeolog primary, HPT/DPT axis
script: /scratch/srussell/srussell/scripts/v9/fig6e_xlaevis_walker_etv1_ucell.R
submitted: $(date -Is)
expected_outputs:
  - /project/def-srussell/nirosha_xeno/figures_v9/fig6/e_panel.pdf
  - /project/def-srussell/nirosha_xeno/figures_v9/fig6/e_panel.png
  - /project/def-srussell/nirosha_xeno/figures_v9/fig6/panel_meta_e.yaml
  - /project/def-srussell/nirosha_xeno/figures_v9/fig6/fig6e_summary.tsv
  - /project/def-srussell/nirosha_xeno/figures_v9/fig6/fig6e_per_sample.tsv
EOF

module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21
export R_LIBS_USER="${HOME}/Rlibs/4.5:${R_LIBS_USER:-}"
cd /project/def-srussell/nirosha_xeno
Rscript /scratch/srussell/srussell/scripts/v9/fig6e_xlaevis_walker_etv1_ucell.R

EXIT=$?
echo "exit_code: ${EXIT}"   >> "$MANIFEST"
echo "finished: $(date -Is)" >> "$MANIFEST"
exit ${EXIT}
