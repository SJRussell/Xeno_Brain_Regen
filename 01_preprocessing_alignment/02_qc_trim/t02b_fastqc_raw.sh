#!/bin/bash
#SBATCH --job-name=t02b_fastqc_raw
#SBATCH --account=def-srussell
#SBATCH --time=01:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=12
#SBATCH --output=/scratch/srussell/srussell/logs/t02b_fastqc_raw_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/t02b_fastqc_raw_%j.err

set -euo pipefail

module purge
module load StdEnv/2023 fastqc/0.12.1 python/3.11

PROJECT_ROOT=/project/def-srussell/nirosha_xeno
FASTQ_DIR="${PROJECT_ROOT}/raw_extracted/fastqs"
QC_RAW_DIR="${PROJECT_ROOT}/qc/raw"
MANIFEST_DIR="${PROJECT_ROOT}/scorial/01_state/slurm_manifests"
mkdir -p "$QC_RAW_DIR" "$MANIFEST_DIR"

cd "$FASTQ_DIR"

# Run FastQC on all merged FASTQs (excluding Undetermined)
echo "[$(date -Iseconds)] Running FastQC on raw FASTQs (parallel via -t)"
FASTQS=$(ls *.fastq.gz | grep -v Undetermined)
N=$(echo "$FASTQS" | wc -l)
echo "  Sample count: $N"

fastqc \
    -o "$QC_RAW_DIR" \
    -t ${SLURM_CPUS_PER_TASK} \
    --noextract \
    --quiet \
    $FASTQS

# Install MultiQC via pip --user if not on cluster modules
if ! command -v multiqc >/dev/null 2>&1; then
    echo "[$(date -Iseconds)] Installing multiqc via pip --user"
    pip install --user multiqc 2>&1 | tail -3
    export PATH="$HOME/.local/bin:$PATH"
fi

# Run MultiQC
echo "[$(date -Iseconds)] Running MultiQC on raw"
cd "$QC_RAW_DIR"
multiqc . -n multiqc_raw -o "${PROJECT_ROOT}/qc" --quiet

# Verify outputs
N_REPORTS=$(ls "$QC_RAW_DIR"/*_fastqc.html 2>/dev/null | wc -l)
echo "[$(date -Iseconds)] FastQC reports generated: $N_REPORTS"
ls -la "${PROJECT_ROOT}/qc/multiqc_raw"* 2>&1 | head -10

EXIT_CODE=$?
MAN="${MANIFEST_DIR}/job_${SLURM_JOB_ID}.yaml"
{
    echo "---"
    echo "job_id: ${SLURM_JOB_ID}"
    echo "job_name: ${SLURM_JOB_NAME}"
    echo "task_id: T-02b"
    echo "stage: stage_A_fastqc_raw_multiqc"
    echo "exit_code: ${EXIT_CODE}"
    echo "end_time: $(date -Iseconds)"
    echo "node: $(hostname)"
    echo "tool: 'fastqc/0.12.1 + multiqc'"
    echo "n_samples: ${N}"
    echo "n_reports: ${N_REPORTS}"
    echo "outputs:"
    echo "  fastqc_raw_dir: ${QC_RAW_DIR}"
    echo "  multiqc_raw: ${PROJECT_ROOT}/qc/multiqc_raw.html"
} > "$MAN"
cat "$MAN"

echo "[$(date -Iseconds)] T-02b complete"
