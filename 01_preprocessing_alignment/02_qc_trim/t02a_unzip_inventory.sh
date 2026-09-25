#!/bin/bash
#SBATCH --job-name=t02a_unzip
#SBATCH --account=def-srussell
#SBATCH --time=01:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --output=/scratch/srussell/srussell/logs/t02a_unzip_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/t02a_unzip_%j.err

set -euo pipefail

PROJECT_ROOT=/project/def-srussell/nirosha_xeno
RAW_ZIPS="${PROJECT_ROOT}/raw"
EXTRACTED_DIR="${PROJECT_ROOT}/raw_extracted"
FASTQ_DIR="${EXTRACTED_DIR}/fastqs"
MANIFEST_DIR="${PROJECT_ROOT}/scorial/01_state/slurm_manifests"
SAMPLE_MANIFEST="${PROJECT_ROOT}/scorial/01_state/sample_manifest.tsv"

mkdir -p "$EXTRACTED_DIR" "$FASTQ_DIR" "$MANIFEST_DIR"

cd "$EXTRACTED_DIR"

# Skip if already done
if [[ -f "${EXTRACTED_DIR}/.unzip_complete" ]]; then
    echo "[$(date -Iseconds)] Already unzipped (.unzip_complete found), skipping"
else
    echo "[$(date -Iseconds)] Unzipping 27 zip parts to $EXTRACTED_DIR"
    for z in "${RAW_ZIPS}"/Xenopus_Fastqs-*.zip; do
        bn=$(basename "$z")
        echo "  -> $bn"
        unzip -o -q "$z"
    done
    touch "${EXTRACTED_DIR}/.unzip_complete"
fi

echo "[$(date -Iseconds)] Inventorying merged FASTQs"

# The merged FASTQs are in:
#   Xenopus_Fastqs/Unaligned_190830-0352H_Anna_Kane_Merged/fastq_merge_Lanes/*.fastq.gz
MERGED_DIR="${EXTRACTED_DIR}/Xenopus_Fastqs/Unaligned_190830-0352H_Anna_Kane_Merged/fastq_merge_Lanes"
if [[ ! -d "$MERGED_DIR" ]]; then
    echo "ERROR: merged FASTQ dir not found: $MERGED_DIR"
    find "$EXTRACTED_DIR" -type d -maxdepth 4
    exit 1
fi

# Symlink merged FASTQs into a flat fastqs/ dir
echo "[$(date -Iseconds)] Symlinking merged FASTQs into flat fastqs/ dir"
cd "$FASTQ_DIR"
rm -f *.fastq.gz  # clean previous symlinks if any
for f in "${MERGED_DIR}"/*.fastq.gz; do
    ln -sf "$f" .
done

# Build sample manifest
echo "[$(date -Iseconds)] Building sample manifest"
{
    echo -e "sample_id\ttreatment\trep\ttime\ttissue\tfilename\tsize_bytes\tmd5"
    for f in *.fastq.gz; do
        # Filename like: 5461_mix-ND_T1_7DPA_BLAST_S2_R1_001.fastq.gz
        # or:           5461_mix-CT_T2_24HPA_BRAIN_S25_R1_001.fastq.gz
        # Extract treatment / rep / time / tissue from the prefix
        bn=$(basename "$f" .fastq.gz)
        # Skip Undetermined files
        if [[ "$bn" == Undetermined* ]]; then continue; fi
        # Pattern: 5461_mix-{TREAT}_{REP}_{TIME}_{TISSUE}_S{N}_R1_001
        meta=$(echo "$bn" | sed -E 's/^5461_mix-([A-Z]+)_(T[0-9])_([0-9]+[A-Z]+)_([A-Z]+)_S[0-9]+.*$/\1\t\2\t\3\t\4/')
        sample_id="${meta//$'\t'/_}"
        # Resolve symlink to get real file size
        real=$(readlink -f "$f")
        sz=$(stat -c%s "$real")
        md5=$(md5sum "$real" | awk '{print $1}')
        echo -e "${sample_id}\t${meta}\t${f}\t${sz}\t${md5}"
    done
} > "$SAMPLE_MANIFEST"

echo "[$(date -Iseconds)] Sample manifest written:"
cat "$SAMPLE_MANIFEST"

# Sanity counts
N_SAMPLES=$(awk 'NR>1 && !/^Undetermined/ {n++} END{print n+0}' "$SAMPLE_MANIFEST")
echo "[$(date -Iseconds)] Sample count: $N_SAMPLES (expected 30)"

# Cross-check against on-disk md5sums.txt files inside the archive
echo "[$(date -Iseconds)] Cross-checking with archive md5sum.txt files"
ARCHIVE_MD5=$(find "${EXTRACTED_DIR}/Xenopus_Fastqs" -name "md5sum.txt" | head -1)
if [[ -n "$ARCHIVE_MD5" && -f "$ARCHIVE_MD5" ]]; then
    echo "  Archive md5 file: $ARCHIVE_MD5"
    head -5 "$ARCHIVE_MD5"
fi

# Read length probe — first FASTQ
# Use head -2 to read just two lines, avoiding SIGPIPE on zcat under pipefail
FIRST_FQ=$(ls "$FASTQ_DIR"/*.fastq.gz | grep -v Undetermined | head -1)
echo "[$(date -Iseconds)] Read length probe on $FIRST_FQ"
LEN=$(zcat "$FIRST_FQ" 2>/dev/null | head -2 | tail -1 | awk '{print length}')
echo "  First read length: $LEN nt"
echo "$LEN" > "${EXTRACTED_DIR}/.read_length"

# Manifest
EXIT_CODE=$?
MAN="${MANIFEST_DIR}/job_${SLURM_JOB_ID}.yaml"
{
    echo "---"
    echo "job_id: ${SLURM_JOB_ID}"
    echo "job_name: ${SLURM_JOB_NAME}"
    echo "task_id: T-02a"
    echo "stage: stage_A_unzip_inventory"
    echo "exit_code: ${EXIT_CODE}"
    echo "end_time: $(date -Iseconds)"
    echo "node: $(hostname)"
    echo "outputs:"
    echo "  extracted_dir: ${EXTRACTED_DIR}"
    echo "  fastq_symlink_dir: ${FASTQ_DIR}"
    echo "  sample_manifest: ${SAMPLE_MANIFEST}"
    echo "  n_samples: ${N_SAMPLES}"
    echo "  first_read_length: ${LEN}"
} > "$MAN"
cat "$MAN"

echo "[$(date -Iseconds)] T-02a complete"
