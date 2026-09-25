#!/bin/bash
#SBATCH --job-name=t01a_dl_xenla
#SBATCH --account=def-srussell
#SBATCH --time=04:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4
#SBATCH --output=/scratch/srussell/srussell/logs/t01a_dl_xenla_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/t01a_dl_xenla_%j.err

set -euo pipefail

PROJECT_ROOT=/project/def-srussell/nirosha_xeno
REF_DIR="${PROJECT_ROOT}/ref/XENLA_10.1"
MANIFEST_DIR="${PROJECT_ROOT}/scorial/01_state/slurm_manifests"
mkdir -p "$REF_DIR" "$MANIFEST_DIR"

cd "$REF_DIR"

XEN_BASE=https://ftp.xenbase.org/pub/Genomics/JGI/Xenla10.1

echo "[$(date -Iseconds)] Starting downloads from $XEN_BASE"

# Files to download: name on Xenbase -> local name
declare -A FILES=(
    ["XENLA_10.1_genome.fa.gz"]="genome.fa.gz"
    ["XENLA_10.1_GCF.gff3.gz"]="annotation_GCF.gff3.gz"
    ["XENLA_10.1_Xenbase.gff3.gz"]="annotation_Xenbase.gff3.gz"
    ["sequences/XENLA_10.1_Xenbase.transcripts.fa.gz"]="transcripts_Xenbase.fa.gz"
)

for src in "${!FILES[@]}"; do
    dst="${FILES[$src]}"
    if [[ -f "$dst" && -s "$dst" ]]; then
        echo "[$(date -Iseconds)] Skipping $dst (already present, $(stat -c%s "$dst") bytes)"
        continue
    fi
    echo "[$(date -Iseconds)] Downloading $src -> $dst"
    curl -fL --retry 3 --retry-delay 10 -o "$dst" "${XEN_BASE}/${src}"
    echo "[$(date -Iseconds)] Done: $dst ($(stat -c%s "$dst") bytes)"
done

# Verify gzip integrity, then decompress (keep originals so md5 can be reverified)
for f in genome.fa.gz annotation_GCF.gff3.gz annotation_Xenbase.gff3.gz transcripts_Xenbase.fa.gz; do
    echo "[$(date -Iseconds)] Verifying gzip integrity of $f"
    gzip -tv "$f"
done

# Decompress (keep .gz copies)
echo "[$(date -Iseconds)] Decompressing genome.fa.gz"
gunzip -kf genome.fa.gz
echo "[$(date -Iseconds)] Decompressing annotation_GCF.gff3.gz"
gunzip -kf annotation_GCF.gff3.gz
echo "[$(date -Iseconds)] Decompressing annotation_Xenbase.gff3.gz"
gunzip -kf annotation_Xenbase.gff3.gz
echo "[$(date -Iseconds)] Decompressing transcripts_Xenbase.fa.gz"
gunzip -kf transcripts_Xenbase.fa.gz

# Generate md5 manifest
echo "[$(date -Iseconds)] Computing md5 sums"
md5sum genome.fa.gz annotation_GCF.gff3.gz annotation_Xenbase.gff3.gz transcripts_Xenbase.fa.gz \
       genome.fa annotation_GCF.gff3 annotation_Xenbase.gff3 transcripts_Xenbase.fa \
       > md5sums.txt
cat md5sums.txt

# Build Salmon decoy file (chromosome names from genome.fa)
echo "[$(date -Iseconds)] Building Salmon decoys.txt"
grep '^>' genome.fa | sed 's/^>//; s/[ \t].*//' > decoys.txt
wc -l decoys.txt

# Concat transcripts + genome for Salmon decoy-aware index
echo "[$(date -Iseconds)] Building gentrome.fa.gz for Salmon"
cat transcripts_Xenbase.fa genome.fa | gzip -1 > gentrome.fa.gz
ls -la gentrome.fa.gz

# Manifest
EXIT_CODE=$?
MAN="${MANIFEST_DIR}/job_${SLURM_JOB_ID}.yaml"
{
    echo "---"
    echo "job_id: ${SLURM_JOB_ID}"
    echo "job_name: ${SLURM_JOB_NAME}"
    echo "task_id: T-01a"
    echo "stage: stage_ref_xenla_10.1_download"
    echo "exit_code: ${EXIT_CODE}"
    echo "end_time: $(date -Iseconds)"
    echo "node: $(hostname)"
    echo "output_files:"
    for f in genome.fa annotation_GCF.gff3 annotation_Xenbase.gff3 transcripts_Xenbase.fa decoys.txt gentrome.fa.gz md5sums.txt; do
        if [[ -f "${REF_DIR}/${f}" ]]; then
            echo "  - path: ${REF_DIR}/${f}"
            echo "    size_bytes: $(stat -c%s "${REF_DIR}/${f}")"
        fi
    done
} > "$MAN"
cat "$MAN"

echo "[$(date -Iseconds)] T-01a complete"
