#!/bin/bash
#SBATCH --job-name=hB_descalzi
#SBATCH --account=def-srussell
#SBATCH --time=01:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=/scratch/srussell/srussell/logs/hB_descalzi_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/hB_descalzi_%j.err

set -euo pipefail
module purge
module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21
export R_LIBS_USER=/home/srussell/Rlibs/4.5

cd /project/def-srussell/nirosha_xeno
Rscript /scratch/srussell/srussell/scripts/hB_descalzi.R

MANIFEST="/project/def-srussell/nirosha_xeno/scorial/01_state/slurm_manifests/job_${SLURM_JOB_ID}.yaml"
cat > "$MANIFEST" <<EOF
---
job_id: ${SLURM_JOB_ID}
job_name: hB_descalzi
task_id: T-09b
stage: comparator_descalzi
purpose: "Descalzi 2017 GSE91396 mouse SNI bulk RNA-seq comparator. limma+trend on log2(RPKM+1) per region (mPFC, NAc, PAG), Hallmark GSEA, gene-level + Hallmark-NES concordance vs Xenopus 7DPA CT brain."
end_time: \$(date -Iseconds)
node: \$(hostname)
exit_code: 0
outputs:
  workspace: /scratch/srussell/srussell/nirosha_xeno_comparators/descalzi
  de_per_region: /project/def-srussell/nirosha_xeno/comparators/descalzi/de_*.tsv
  gsea_per_region: /project/def-srussell/nirosha_xeno/comparators/descalzi/gsea_*.tsv
  concordance: /project/def-srussell/nirosha_xeno/comparators/descalzi/xenopus_concordance.tsv
EOF
echo "Manifest written: $MANIFEST"
