#!/bin/bash
#SBATCH --job-name=t04_stage_C
#SBATCH --account=def-srussell
#SBATCH --time=01:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --output=/scratch/srussell/srussell/logs/t04_stage_C_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/t04_stage_C_%j.err

set -euo pipefail

module purge
module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21

PROJECT_ROOT=/project/def-srussell/nirosha_xeno
MANIFEST_DIR="${PROJECT_ROOT}/scorial/01_state/slurm_manifests"
mkdir -p "$MANIFEST_DIR"

cd "${PROJECT_ROOT}"

echo "[$(date -Iseconds)] Running Stage C R script"
echo "  R version: $(R --version | head -1)"

Rscript /scratch/srussell/srussell/scripts/t04_stage_C_qc.R

EXIT_CODE=$?
echo "[$(date -Iseconds)] Stage C exited ${EXIT_CODE}"

MAN="${MANIFEST_DIR}/job_${SLURM_JOB_ID}.yaml"
{
    echo "---"
    echo "job_id: ${SLURM_JOB_ID}"
    echo "job_name: ${SLURM_JOB_NAME}"
    echo "task_id: T-04"
    echo "stage: stage_C_qc_outlier_triage"
    echo "exit_code: ${EXIT_CODE}"
    echo "end_time: $(date -Iseconds)"
    echo "node: $(hostname)"
    echo "tool: 'R/4.5.0 + r-bundle-bioconductor/3.21 (edgeR + DESeq2 + limma)'"
    echo "outputs:"
    echo "  figures: ${PROJECT_ROOT}/figures/stage_C/"
    echo "  data:    ${PROJECT_ROOT}/deg/"
    echo "    - sample_metadata.{tsv,rds}"
    echo "    - dgelist_filtered_tmm.rds"
    echo "    - vsd_blind.rds"
    echo "    - pca_BRAIN.tsv, pca_BLAST.tsv"
    echo "    - outlier_diagnostics.tsv"
} > "$MAN"
cat "$MAN"

echo "[$(date -Iseconds)] T-04 complete"
