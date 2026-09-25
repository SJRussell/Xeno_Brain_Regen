#!/bin/bash
#SBATCH --job-name=t02d2_fastqc_trim
#SBATCH --account=def-srussell
#SBATCH --time=00:45:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=12
#SBATCH --output=/scratch/srussell/srussell/logs/t02d2_fastqc_trim_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/t02d2_fastqc_trim_%j.err

set -euo pipefail

module purge
module load StdEnv/2023 fastqc/0.12.1

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

N_REPORTS=$(ls "$QC_TRIM_DIR"/*_fastqc.html 2>/dev/null | wc -l)
echo "[$(date -Iseconds)] FastQC reports generated: $N_REPORTS"

EXIT_CODE=$?
MAN="${MANIFEST_DIR}/job_${SLURM_JOB_ID}.yaml"
{
    echo "---"
    echo "job_id: ${SLURM_JOB_ID}"
    echo "job_name: ${SLURM_JOB_NAME}"
    echo "task_id: T-02d"
    echo "stage: stage_A_fastqc_trim"
    echo "exit_code: ${EXIT_CODE}"
    echo "end_time: $(date -Iseconds)"
    echo "node: $(hostname)"
    echo "tool: fastqc/0.12.1"
    echo "n_samples: ${N}"
    echo "n_reports: ${N_REPORTS}"
    echo "outputs:"
    echo "  fastqc_trim_dir: ${QC_TRIM_DIR}"
    echo "note: 'multiqc aggregation runs separately on login node via ~/.venvs/multiqc'"
} > "$MAN"
cat "$MAN"

echo "[$(date -Iseconds)] T-02d FastQC complete (multiqc to follow on login node)"
