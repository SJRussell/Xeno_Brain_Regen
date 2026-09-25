#!/bin/bash
#SBATCH --job-name=t02d_fastqc_trim
#SBATCH --account=def-srussell
#SBATCH --time=01:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=12
#SBATCH --output=/scratch/srussell/srussell/logs/t02d_fastqc_trim_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/t02d_fastqc_trim_%j.err

set -euo pipefail

module purge
module load StdEnv/2023 fastqc/0.12.1 python/3.11

PROJECT_ROOT=/project/def-srussell/nirosha_xeno
TRIM_DIR="${PROJECT_ROOT}/trim"
QC_TRIM_DIR="${PROJECT_ROOT}/qc/trim"
MANIFEST_DIR="${PROJECT_ROOT}/scorial/01_state/slurm_manifests"
mkdir -p "$QC_TRIM_DIR" "$MANIFEST_DIR"

cd "$TRIM_DIR"

echo "[$(date -Iseconds)] Running FastQC on trimmed FASTQs"
FASTQS=$(ls *.trim.fastq.gz)
N=$(echo "$FASTQS" | wc -l)
echo "  Trimmed sample count: $N"

fastqc \
    -o "$QC_TRIM_DIR" \
    -t ${SLURM_CPUS_PER_TASK} \
    --noextract \
    --quiet \
    $FASTQS

# Ensure multiqc is on PATH (installed in T-02b)
export PATH="$HOME/.local/bin:$PATH"

# MultiQC across both raw + trim + fastp
echo "[$(date -Iseconds)] Running MultiQC across qc/raw, qc/trim, qc/fastp"
cd "${PROJECT_ROOT}/qc"
multiqc raw trim fastp -n multiqc_combined -o "${PROJECT_ROOT}/qc" --force --quiet

multiqc trim fastp -n multiqc_trim -o "${PROJECT_ROOT}/qc" --force --quiet

N_REPORTS=$(ls "$QC_TRIM_DIR"/*_fastqc.html 2>/dev/null | wc -l)

EXIT_CODE=$?
MAN="${MANIFEST_DIR}/job_${SLURM_JOB_ID}.yaml"
{
    echo "---"
    echo "job_id: ${SLURM_JOB_ID}"
    echo "job_name: ${SLURM_JOB_NAME}"
    echo "task_id: T-02d"
    echo "stage: stage_A_fastqc_trim_multiqc"
    echo "exit_code: ${EXIT_CODE}"
    echo "end_time: $(date -Iseconds)"
    echo "node: $(hostname)"
    echo "tool: 'fastqc/0.12.1 + multiqc'"
    echo "n_samples: ${N}"
    echo "n_reports: ${N_REPORTS}"
    echo "outputs:"
    echo "  fastqc_trim_dir: ${QC_TRIM_DIR}"
    echo "  multiqc_trim: ${PROJECT_ROOT}/qc/multiqc_trim.html"
    echo "  multiqc_combined: ${PROJECT_ROOT}/qc/multiqc_combined.html"
} > "$MAN"
cat "$MAN"

echo "[$(date -Iseconds)] T-02d complete (Stage A done)"
