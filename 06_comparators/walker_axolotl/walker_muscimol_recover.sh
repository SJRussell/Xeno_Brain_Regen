#!/bin/bash
#SBATCH --job-name=walker_musc_recover
#SBATCH --account=def-srussell
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=/scratch/srussell/srussell/logs/walker_musc_recover_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/walker_musc_recover_%j.err

# RED-1 (per Stewart's iter-19 correction of DEC-019).
# Re-derive SRR31800874 (Walker muscimol arm) FASTQs via NCBI prefetch + fasterq-dump.
# The previous ENA-mirror extract was reproducibly corrupt at 204M barcodes
# (alevin failed twice on 13133998_3, 13136766_3 — gzip -t passed, FASTQ-record
# corruption). prefetch pulls the .sra archive directly; fasterq-dump rebuilds
# FASTQs locally — different derivation from the ENA mirror.

set -euo pipefail

module purge
module load StdEnv/2023 sra-toolkit/3.0.9

SRR="SRR31800874"
STAGE="/scratch/srussell/srussell/nirosha_xeno_comparators/walker/raw/muscimol_brain_recover"
DEST="/scratch/srussell/srussell/nirosha_xeno_comparators/walker/raw/muscimol_brain"
mkdir -p "$STAGE"
cd "$STAGE"

echo "=========================================================="
echo "[$(date -Iseconds)] RED-1 muscimol recovery for $SRR"
echo "  Stage: $STAGE"
echo "  Final dest: $DEST"
echo "  Old (corrupt) R2: $DEST/${SRR}_2.fastq.gz (will be preserved as .ena_corrupt)"
echo "=========================================================="

# Step 1: prefetch — pulls .sra archive from NCBI (different source than ENA mirror)
echo "[$(date -Iseconds)] Step 1: prefetch $SRR ..."
# Allow up to ~50 GB for the .sra archive
prefetch --max-size 100g --output-directory "$STAGE" "$SRR"
ls -lh "$STAGE/$SRR/"
SRA_FILE="$STAGE/$SRR/${SRR}.sra"
[ -f "$SRA_FILE" ] || { echo "ERROR: prefetch did not produce expected .sra"; exit 2; }
echo "[$(date -Iseconds)] prefetch complete; .sra size: $(stat -c %s "$SRA_FILE")"

# Step 2: fasterq-dump --split-files
echo "[$(date -Iseconds)] Step 2: fasterq-dump --split-files ..."
TMPSPACE="$STAGE/_fasterq_tmp"
mkdir -p "$TMPSPACE"
fasterq-dump \
    --split-files \
    --threads "$SLURM_CPUS_PER_TASK" \
    --temp "$TMPSPACE" \
    --outdir "$STAGE" \
    "$SRA_FILE"
rm -rf "$TMPSPACE"

R1_RAW="$STAGE/${SRR}_1.fastq"
R2_RAW="$STAGE/${SRR}_2.fastq"
[ -f "$R1_RAW" ] || { echo "ERROR: missing R1 from fasterq-dump"; exit 3; }
[ -f "$R2_RAW" ] || { echo "ERROR: missing R2 from fasterq-dump"; exit 3; }
echo "[$(date -Iseconds)] fasterq-dump complete:"
ls -lh "$R1_RAW" "$R2_RAW"

# Step 3: gzip (use pigz if available, else gzip)
echo "[$(date -Iseconds)] Step 3: compress with pigz ..."
if command -v pigz >/dev/null 2>&1; then
    pigz -p "$SLURM_CPUS_PER_TASK" "$R1_RAW" "$R2_RAW"
else
    gzip "$R1_RAW" "$R2_RAW"
fi

R1_GZ="${R1_RAW}.gz"
R2_GZ="${R2_RAW}.gz"
ls -lh "$R1_GZ" "$R2_GZ"

# Step 4: gzip integrity test
echo "[$(date -Iseconds)] Step 4: gzip -t validation ..."
gzip -t "$R1_GZ" && echo "  R1 OK" || { echo "  R1 FAIL"; exit 4; }
gzip -t "$R2_GZ" && echo "  R2 OK" || { echo "  R2 FAIL"; exit 4; }

# Step 5: compare to old corrupt R2 — if MD5 matches, re-derivation produced the same file
OLD_R2="$DEST/${SRR}_2.fastq.gz"
if [ -f "$OLD_R2" ]; then
    echo "[$(date -Iseconds)] Step 5: comparing new R2 to old (ENA-mirror) R2 ..."
    NEW_SIZE=$(stat -c %s "$R2_GZ")
    OLD_SIZE=$(stat -c %s "$OLD_R2")
    echo "  Old R2 size: $OLD_SIZE"
    echo "  New R2 size: $NEW_SIZE"
    NEW_MD5=$(md5sum "$R2_GZ" | awk '{print $1}')
    OLD_MD5=$(md5sum "$OLD_R2" | awk '{print $1}')
    echo "  Old R2 MD5:  $OLD_MD5"
    echo "  New R2 MD5:  $NEW_MD5"
    if [ "$NEW_MD5" = "$OLD_MD5" ]; then
        echo "  WARNING: identical bytes — recovery via prefetch produced the same content."
        echo "  This means corruption is in the .sra archive itself, not the ENA mirror."
        echo "  Will need AWS S3 sra-src bucket fallback."
    else
        echo "  Files differ — prefetch produced an alternate derivation."
    fi
fi

# Step 6: stage the new files into DEST as the canonical FASTQs
# Preserve the old (corrupt) R2 alongside.
echo "[$(date -Iseconds)] Step 6: staging into $DEST ..."
mkdir -p "$DEST"
if [ -f "$DEST/${SRR}_2.fastq.gz" ]; then
    mv -v "$DEST/${SRR}_2.fastq.gz" "$DEST/${SRR}_2.fastq.gz.ena_corrupt"
fi
if [ -f "$DEST/${SRR}_1.fastq.gz" ]; then
    # R1 was not corrupt before; back it up only if size differs
    OLD_R1_SIZE=$(stat -c %s "$DEST/${SRR}_1.fastq.gz")
    NEW_R1_SIZE=$(stat -c %s "$R1_GZ")
    if [ "$OLD_R1_SIZE" != "$NEW_R1_SIZE" ]; then
        mv -v "$DEST/${SRR}_1.fastq.gz" "$DEST/${SRR}_1.fastq.gz.ena_orig"
    else
        rm -v "$DEST/${SRR}_1.fastq.gz"
    fi
fi
mv -v "$R1_GZ" "$DEST/"
mv -v "$R2_GZ" "$DEST/"

echo "[$(date -Iseconds)] Final DEST contents:"
ls -lh "$DEST/"

# Manifest
PROJ="/project/def-srussell/nirosha_xeno"
MANIFEST="$PROJ/scorial/01_state/slurm_manifests/job_${SLURM_JOB_ID}.yaml"
cat > "$MANIFEST" <<EOF
---
job_id: ${SLURM_JOB_ID}
job_name: walker_musc_recover
task_id: RED-1
stage: comparator_walker_muscimol_recovery
purpose: "RED-1 (iter-19): re-derive SRR31800874 FASTQs via NCBI prefetch (different source than ENA mirror) to recover muscimol arm dropped in DEC-019."
end_time: $(date -Iseconds)
node: $(hostname)
exit_code: 0
related_decisions: [DEC-019]
outputs:
  - $DEST/${SRR}_1.fastq.gz
  - $DEST/${SRR}_2.fastq.gz
old_files_preserved:
  - $DEST/${SRR}_2.fastq.gz.ena_corrupt
EOF
echo "Manifest: $MANIFEST"
echo "[$(date -Iseconds)] RED-1 muscimol recovery complete."
