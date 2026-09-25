#!/bin/bash
#SBATCH --job-name=t05_stage_D
#SBATCH --account=def-srussell
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=/scratch/srussell/srussell/logs/t05_stage_D_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/t05_stage_D_%j.err

set -euo pipefail

module purge
module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21

PROJECT_ROOT=/project/def-srussell/nirosha_xeno
MANIFEST_DIR="${PROJECT_ROOT}/scorial/01_state/slurm_manifests"
mkdir -p "$MANIFEST_DIR"

cd "${PROJECT_ROOT}"

echo "[$(date -Iseconds)] Running Stage D R script"
Rscript /scratch/srussell/srussell/scripts/t05_stage_D_de.R

EXIT_CODE=$?
echo "[$(date -Iseconds)] Stage D exited ${EXIT_CODE}"

MAN="${MANIFEST_DIR}/job_${SLURM_JOB_ID}.yaml"
{
    echo "---"
    echo "job_id: ${SLURM_JOB_ID}"
    echo "job_name: ${SLURM_JOB_NAME}"
    echo "task_id: T-05"
    echo "stage: stage_D_de"
    echo "exit_code: ${EXIT_CODE}"
    echo "end_time: $(date -Iseconds)"
    echo "node: $(hostname)"
    echo "tool: 'limma + voomWithQualityWeights (primary), DESeq2 LRT (interaction confirmatory), edgeR'"
    echo "outputs:"
    echo "  - deg/voom_BRAIN_primary.{rds,xlsx}"
    echo "  - deg/voom_BRAIN_sensitivity.{rds,xlsx}"
    echo "  - deg/voom_BLAST_primary.{rds,xlsx}"
    echo "  - deg/deseq2_lrt_*.rds"
    echo "  - deg/all_de_summary.tsv"
} > "$MAN"
cat "$MAN"
