#!/bin/bash
#SBATCH --job-name=t02c_fastp
#SBATCH --account=def-srussell
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --output=/scratch/srussell/srussell/logs/t02c_fastp_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/t02c_fastp_%j.err

set -euo pipefail

module purge
module load StdEnv/2023 fastp/0.24.0

PROJECT_ROOT=/project/def-srussell/nirosha_xeno
FASTQ_DIR="${PROJECT_ROOT}/raw_extracted/fastqs"
TRIM_DIR="${PROJECT_ROOT}/trim"
TRIM_REPORTS="${PROJECT_ROOT}/qc/fastp"
MANIFEST_DIR="${PROJECT_ROOT}/scorial/01_state/slurm_manifests"
mkdir -p "$TRIM_DIR" "$TRIM_REPORTS" "$MANIFEST_DIR"

cd "$FASTQ_DIR"

# Run fastp single-end with default Illumina adapters, length>=36, Q>=20.
# Run 4 samples concurrently, each using 4 threads (= 16 total CPUs).

run_fastp() {
    local in_fq="$1"
    local bn=$(basename "$in_fq" .fastq.gz)
    local out_fq="${TRIM_DIR}/${bn}.trim.fastq.gz"
    local html="${TRIM_REPORTS}/${bn}.fastp.html"
    local json="${TRIM_REPORTS}/${bn}.fastp.json"

    if [[ -f "$out_fq" && -s "$out_fq" ]]; then
        echo "  [SKIP] $bn (already trimmed)"
        return 0
    fi

    fastp \
        --in1 "$in_fq" \
        --out1 "$out_fq" \
        --html "$html" \
        --json "$json" \
        --length_required 36 \
        --qualified_quality_phred 20 \
        --thread 4 \
        2>&1 | tail -10
    echo "  [DONE] $bn"
}
export -f run_fastp
export TRIM_DIR TRIM_REPORTS

FASTQS=$(ls *.fastq.gz | grep -v Undetermined)
N=$(echo "$FASTQS" | wc -l)
echo "[$(date -Iseconds)] Running fastp on $N samples (4-way parallel, 4 threads each)"

# Run 4 in parallel using xargs
echo "$FASTQS" | xargs -n 1 -P 4 -I {} bash -c 'run_fastp "$@"' _ {}

# Verify
N_TRIMMED=$(ls "$TRIM_DIR"/*.trim.fastq.gz 2>/dev/null | wc -l)
N_JSON=$(ls "$TRIM_REPORTS"/*.fastp.json 2>/dev/null | wc -l)

echo "[$(date -Iseconds)] Trimmed FASTQs: $N_TRIMMED"
echo "[$(date -Iseconds)] fastp JSON reports: $N_JSON"

EXIT_CODE=$?
MAN="${MANIFEST_DIR}/job_${SLURM_JOB_ID}.yaml"
{
    echo "---"
    echo "job_id: ${SLURM_JOB_ID}"
    echo "job_name: ${SLURM_JOB_NAME}"
    echo "task_id: T-02c"
    echo "stage: stage_A_fastp_trim"
    echo "exit_code: ${EXIT_CODE}"
    echo "end_time: $(date -Iseconds)"
    echo "node: $(hostname)"
    echo "tool: fastp/0.24.0"
    echo "parameters:"
    echo "  length_required: 36"
    echo "  qualified_quality_phred: 20"
    echo "  adapters: default Illumina (auto-detect)"
    echo "  parallel_samples: 4"
    echo "  threads_per_sample: 4"
    echo "n_input: ${N}"
    echo "n_trimmed: ${N_TRIMMED}"
    echo "outputs:"
    echo "  trim_dir: ${TRIM_DIR}"
    echo "  fastp_reports_dir: ${TRIM_REPORTS}"
} > "$MAN"
cat "$MAN"

echo "[$(date -Iseconds)] T-02c complete"
