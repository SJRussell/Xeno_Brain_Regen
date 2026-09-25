#!/bin/bash
#SBATCH --job-name=t01c_salmon_idx
#SBATCH --account=def-srussell
#SBATCH --time=02:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=8
#SBATCH --output=/scratch/srussell/srussell/logs/t01c_salmon_idx_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/t01c_salmon_idx_%j.err

set -euo pipefail

module purge
module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 salmon/1.10.2

PROJECT_ROOT=/project/def-srussell/nirosha_xeno
REF_DIR="${PROJECT_ROOT}/ref/XENLA_10.1"
INDEX_DIR="${REF_DIR}/salmon_index"
MANIFEST_DIR="${PROJECT_ROOT}/scorial/01_state/slurm_manifests"
mkdir -p "$INDEX_DIR" "$MANIFEST_DIR"

cd "$REF_DIR"

# Decoy-aware Salmon index. k=25 matches pipeline_plan default.
# gentrome.fa.gz (transcripts + genome) and decoys.txt built in T-01a.
echo "[$(date -Iseconds)] Starting Salmon index"
echo "  Gentrome: ${REF_DIR}/gentrome.fa.gz"
echo "  Decoys:   ${REF_DIR}/decoys.txt ($(wc -l < ${REF_DIR}/decoys.txt) entries)"

salmon index \
    -t "${REF_DIR}/gentrome.fa.gz" \
    -d "${REF_DIR}/decoys.txt" \
    -i "${INDEX_DIR}" \
    -k 25 \
    -p ${SLURM_CPUS_PER_TASK}

EXIT_CODE=$?
echo "[$(date -Iseconds)] Salmon exited ${EXIT_CODE}"

ls -la "${INDEX_DIR}"

MAN="${MANIFEST_DIR}/job_${SLURM_JOB_ID}.yaml"
{
    echo "---"
    echo "job_id: ${SLURM_JOB_ID}"
    echo "job_name: ${SLURM_JOB_NAME}"
    echo "task_id: T-01c"
    echo "stage: stage_ref_xenla_10.1_salmon_index"
    echo "exit_code: ${EXIT_CODE}"
    echo "end_time: $(date -Iseconds)"
    echo "node: $(hostname)"
    echo "tool: salmon/1.10.2"
    echo "parameters:"
    echo "  k: 25"
    echo "  decoy_aware: true"
    echo "  threads: ${SLURM_CPUS_PER_TASK}"
    echo "output_files:"
    if [[ -d "${INDEX_DIR}" ]]; then
        for f in info.json mphf.bin pos.bin seq.bin ref_k25_fixed.bin; do
            if [[ -f "${INDEX_DIR}/${f}" ]]; then
                echo "  - path: ${INDEX_DIR}/${f}"
                echo "    size_bytes: $(stat -c%s "${INDEX_DIR}/${f}")"
            fi
        done
    fi
} > "$MAN"
cat "$MAN"

echo "[$(date -Iseconds)] T-01c complete"
