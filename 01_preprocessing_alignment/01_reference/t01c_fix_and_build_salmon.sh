#!/bin/bash
#SBATCH --job-name=t01c_fix_salmon
#SBATCH --account=def-srussell
#SBATCH --time=02:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=8
#SBATCH --output=/scratch/srussell/srussell/logs/t01c_fix_salmon_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/t01c_fix_salmon_%j.err

set -euo pipefail

module purge
module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 salmon/1.10.2

PROJECT_ROOT=/project/def-srussell/nirosha_xeno
REF_DIR="${PROJECT_ROOT}/ref/XENLA_10.1"
INDEX_DIR="${REF_DIR}/salmon_index"
MANIFEST_DIR="${PROJECT_ROOT}/scorial/01_state/slurm_manifests"
mkdir -p "$INDEX_DIR" "$MANIFEST_DIR"

cd "$REF_DIR"

# Step 1: filter transcripts FASTA to remove records with empty headers
# or sequences shorter than 31 nt (Salmon k=25 plus a margin).
# Xenbase ships ~1124 malformed bare ">" records that break Salmon indexing.
echo "[$(date -Iseconds)] Filtering transcripts FASTA"

awk '
    BEGIN { name=""; seq="" }
    /^>/ {
        if (name != "" && length(seq) >= 31) {
            print name
            print seq
        }
        name = $0
        seq = ""
        # treat header as empty if it is just ">"
        if (name == ">") name = ""
        next
    }
    { seq = seq $0 }
    END {
        if (name != "" && length(seq) >= 31) {
            print name
            print seq
        }
    }
' transcripts_Xenbase.fa > transcripts_Xenbase.clean.fa

ORIG=$(grep -c '^>' transcripts_Xenbase.fa)
CLEAN=$(grep -c '^>' transcripts_Xenbase.clean.fa)
DROPPED=$((ORIG - CLEAN))
echo "  Original transcripts: $ORIG"
echo "  Clean transcripts:    $CLEAN"
echo "  Dropped:              $DROPPED"

# Sanity: ensure no duplicate transcript IDs (first whitespace-delimited token)
DUP=$(grep '^>' transcripts_Xenbase.clean.fa | awk '{print $1}' | sort | uniq -d | head)
if [[ -n "$DUP" ]]; then
    echo "WARN: duplicate transcript IDs found:"
    echo "$DUP" | head -5
fi

# Step 2: rebuild gentrome
echo "[$(date -Iseconds)] Rebuilding gentrome.fa.gz"
cat transcripts_Xenbase.clean.fa genome.fa | gzip -1 > gentrome.fa.gz
ls -la gentrome.fa.gz

# Step 3: build Salmon index
echo "[$(date -Iseconds)] Starting Salmon index"
echo "  Gentrome: ${REF_DIR}/gentrome.fa.gz"
echo "  Decoys:   ${REF_DIR}/decoys.txt ($(wc -l < ${REF_DIR}/decoys.txt) entries)"

# Clean any partial index from previous failed attempt
rm -rf "${INDEX_DIR}"
mkdir -p "${INDEX_DIR}"

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
    echo "preprocessing:"
    echo "  filter: 'awk drop empty-header and seq<31nt'"
    echo "  transcripts_dropped: ${DROPPED}"
    echo "  transcripts_kept: ${CLEAN}"
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
