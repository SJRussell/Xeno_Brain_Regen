#!/bin/bash
#SBATCH --job-name=walker_musc_finish
#SBATCH --account=def-srussell
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=/scratch/srussell/srussell/logs/walker_musc_finish_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/walker_musc_finish_%j.err

# RED-1 finish: pigz R2 (R1 already compressed), then move both to dest.
# Picks up from cancelled walker_musc_recover (job 13149743) which had:
#   - SRR31800874_1.fastq.gz (40.7 GB, complete)
#   - SRR31800874_2.fastq    (182 GB raw, awaiting compression)
#   - SRR31800874_2.fastq.gz (partial — to be deleted)
set -euo pipefail

STAGE="/scratch/srussell/srussell/nirosha_xeno_comparators/walker/raw/muscimol_brain_recover"
DEST="/scratch/srussell/srussell/nirosha_xeno_comparators/walker/raw/muscimol_brain"
SRR="SRR31800874"
cd "$STAGE"

R1_GZ="$STAGE/${SRR}_1.fastq.gz"
R2_RAW="$STAGE/${SRR}_2.fastq"
R2_GZ="$STAGE/${SRR}_2.fastq.gz"

echo "[$(date -Iseconds)] RED-1 finish: pigz R2 + move to dest"
echo "  R1.fastq.gz exists: $(ls -lh "$R1_GZ" 2>&1)"
echo "  R2.fastq exists:    $(ls -lh "$R2_RAW" 2>&1)"
echo "  R2.fastq.gz partial (will be removed): $(ls -lh "$R2_GZ" 2>&1 || echo none)"

# Sanity: R1.fastq.gz must be intact
[ -f "$R1_GZ" ] || { echo "ERROR: R1.fastq.gz missing"; exit 2; }
gzip -t "$R1_GZ" || { echo "ERROR: R1.fastq.gz corrupt"; exit 3; }

# Remove any partial R2 .gz from interrupted pigz
[ -f "$R2_GZ" ] && rm -v "$R2_GZ"

# Verify R2 raw is intact
[ -f "$R2_RAW" ] || { echo "ERROR: R2.fastq raw missing — fasterq-dump output lost?"; exit 4; }

# Compress R2 with pigz
echo "[$(date -Iseconds)] Compressing R2 with pigz..."
pigz -p "$SLURM_CPUS_PER_TASK" "$R2_RAW"
[ -f "$R2_GZ" ] || { echo "ERROR: pigz did not produce R2.fastq.gz"; exit 5; }

echo "[$(date -Iseconds)] pigz complete:"
ls -lh "$R1_GZ" "$R2_GZ"

# Validate
echo "[$(date -Iseconds)] Validating gzip integrity..."
gzip -t "$R1_GZ" && echo "  R1 OK" || { echo "  R1 FAIL"; exit 6; }
gzip -t "$R2_GZ" && echo "  R2 OK" || { echo "  R2 FAIL"; exit 6; }

# Compare to old corrupt R2
OLD_R2="$DEST/${SRR}_2.fastq.gz"
if [ -f "$OLD_R2" ]; then
    echo "[$(date -Iseconds)] Comparing new R2 to old (ENA-mirror) R2..."
    NEW_SIZE=$(stat -c %s "$R2_GZ")
    OLD_SIZE=$(stat -c %s "$OLD_R2")
    NEW_MD5=$(md5sum "$R2_GZ" | awk '{print $1}')
    OLD_MD5=$(md5sum "$OLD_R2" | awk '{print $1}')
    echo "  Old R2: size=$OLD_SIZE  md5=$OLD_MD5"
    echo "  New R2: size=$NEW_SIZE  md5=$NEW_MD5"
    if [ "$NEW_MD5" = "$OLD_MD5" ]; then
        echo "  WARNING: identical — recovery via prefetch produced same content; corruption is in archive not mirror."
    else
        echo "  Files differ — prefetch produced an alternate derivation (expected)."
    fi
fi

# Move into DEST as canonical
echo "[$(date -Iseconds)] Staging into $DEST..."
mkdir -p "$DEST"
if [ -f "$DEST/${SRR}_2.fastq.gz" ]; then
    mv -v "$DEST/${SRR}_2.fastq.gz" "$DEST/${SRR}_2.fastq.gz.ena_corrupt"
fi
if [ -f "$DEST/${SRR}_1.fastq.gz" ]; then
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

PROJ="/project/def-srussell/nirosha_xeno"
MANIFEST="$PROJ/scorial/01_state/slurm_manifests/job_${SLURM_JOB_ID}.yaml"
cat > "$MANIFEST" <<EOF
---
job_id: ${SLURM_JOB_ID}
job_name: walker_musc_finish
task_id: RED-1-finish
stage: comparator_walker_muscimol_finish_after_timeout
purpose: "RED-1 follow-up after walker_musc_recover (13149743) cancelled near 8h SLURM limit. Resumes pigz on R2 + moves canonical fastqs to DEST. R1 was already compressed by the cancelled job."
end_time: $(date -Iseconds)
node: $(hostname)
exit_code: 0
related_decisions: [DEC-019]
related_jobs: [13149743]
outputs:
  - $DEST/${SRR}_1.fastq.gz
  - $DEST/${SRR}_2.fastq.gz
old_files_preserved:
  - $DEST/${SRR}_2.fastq.gz.ena_corrupt
EOF
echo "Manifest: $MANIFEST"
echo "[$(date -Iseconds)] RED-1 finish complete."
