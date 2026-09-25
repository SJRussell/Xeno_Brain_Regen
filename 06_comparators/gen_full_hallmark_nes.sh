#!/bin/bash
#SBATCH --job-name=v9_gen_full_hm
#SBATCH --account=def-srussell
#SBATCH --time=00:20:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=/scratch/srussell/srussell/logs/v9_gen_full_hm_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/v9_gen_full_hm_%j.err

set -euo pipefail
MANIFEST="/project/def-srussell/nirosha_xeno/scorial/01_state/slurm_manifests/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.yaml"
cat > "$MANIFEST" <<EOF
job_id: ${SLURM_JOB_ID}
job_name: ${SLURM_JOB_NAME}
task: v9 Group C fix - generate full 50-Hallmark NES TSVs (X. laevis extract + Descalzi/Erdei fgsea no-filter)
script: /scratch/srussell/srussell/scripts/v9/gen_full_hallmark_nes.R
submitted: $(date -Is)
EOF

module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21
export R_LIBS_USER="${HOME}/Rlibs/4.5:${R_LIBS_USER:-}"
cd /project/def-srussell/nirosha_xeno
Rscript /scratch/srussell/srussell/scripts/v9/gen_full_hallmark_nes.R

EXIT=$?
echo "exit_code: ${EXIT}"   >> "$MANIFEST"
echo "finished: $(date -Is)" >> "$MANIFEST"
exit ${EXIT}
