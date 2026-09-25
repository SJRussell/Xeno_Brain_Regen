#!/bin/bash
#SBATCH --job-name=t08a_mofa
#SBATCH --account=def-srussell
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=/scratch/srussell/srussell/logs/t08a_mofa_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/t08a_mofa_%j.err

set -euo pipefail
module purge
module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21
export R_LIBS_USER=/home/srussell/Rlibs/4.5
export OMP_NUM_THREADS=4

cd /project/def-srussell/nirosha_xeno
Rscript /scratch/srussell/srussell/scripts/t08a_stage_G_mofa.R

# Manifest
MANIFEST="/project/def-srussell/nirosha_xeno/scorial/01_state/slurm_manifests/job_${SLURM_JOB_ID}.yaml"
cat > "$MANIFEST" <<EOF
---
job_id: ${SLURM_JOB_ID}
job_name: t08a_mofa
task_id: T-08a
stage: stage_G_mofa
purpose: "MOFA+ training on n=12 paired BRAIN+BLAST views (24HPA + 7DPA × CT/ND × 3 reps), top 8000 variable genes per view, K=10, ARD+spikeslab."
end_time: $(date -Iseconds)
node: $(hostname)
exit_code: 0
outputs:
  workspace: /scratch/srussell/srussell/nirosha_xeno_stage_G/mofa
  hdf5: /scratch/srussell/srussell/nirosha_xeno_stage_G/mofa/mofa_model.hdf5
  factor_variance_explained: /project/def-srussell/nirosha_xeno/mofa/factor_variance_explained.tsv
  factor_scores: /project/def-srussell/nirosha_xeno/mofa/factor_scores.tsv
  factor_loadings_top: /project/def-srussell/nirosha_xeno/mofa/factor_loadings_top.tsv
  factor_treatment_associations: /project/def-srussell/nirosha_xeno/mofa/factor_treatment_associations.tsv
EOF
echo "Manifest written: $MANIFEST"
