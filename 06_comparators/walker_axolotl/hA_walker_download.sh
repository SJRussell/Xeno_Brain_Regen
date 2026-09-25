#!/bin/bash
#SBATCH --job-name=hA_walker_dl
#SBATCH --account=def-srussell
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8
#SBATCH --output=/scratch/srussell/srussell/logs/hA_walker_dl_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/hA_walker_dl_%j.err

# H-A Walker comparator (DEC-017): download raw 10X scRNA FASTQs from ENA.
# 4 samples × paired R1+R2 .fastq.gz from ftp.sra.ebi.ac.uk:
#   SRR31800877 KE1 Uninjured Brain (~125 GB pair)
#   SRR31800876 KE2 Tail Brain      (~126 GB pair)
#   SRR31800875 Limb Brain          (~112 GB pair)
#   SRR31800874 Muscimol Brain      (~122 GB pair)
# Total ~485 GB. ENA FTP direct is faster than NCBI prefetch+fasterq-dump.

set -euo pipefail
WORKROOT="/scratch/srussell/srussell/nirosha_xeno_comparators/walker"
RAW="$WORKROOT/raw"
mkdir -p "$RAW"
cd "$RAW"

# Direct ENA URLs from filereport API (subdir is "0" + last 2 digits, NOT last 3 chars):
#   SRR31800877 -> /vol1/fastq/SRR318/077/SRR31800877/   (uninjured_brain)
#   SRR31800876 -> /vol1/fastq/SRR318/076/SRR31800876/   (tail_brain)
#   SRR31800875 -> /vol1/fastq/SRR318/075/SRR31800875/   (limb_brain)
#   SRR31800874 -> /vol1/fastq/SRR318/074/SRR31800874/   (muscimol_brain)
URLS_R1=(
  "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR318/077/SRR31800877/SRR31800877_1.fastq.gz"
  "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR318/076/SRR31800876/SRR31800876_1.fastq.gz"
  "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR318/075/SRR31800875/SRR31800875_1.fastq.gz"
  "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR318/074/SRR31800874/SRR31800874_1.fastq.gz"
)
URLS_R2=(
  "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR318/077/SRR31800877/SRR31800877_2.fastq.gz"
  "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR318/076/SRR31800876/SRR31800876_2.fastq.gz"
  "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR318/075/SRR31800875/SRR31800875_2.fastq.gz"
  "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR318/074/SRR31800874/SRR31800874_2.fastq.gz"
)
LABELS=(uninjured_brain tail_brain limb_brain muscimol_brain)
SRRS=(SRR31800877 SRR31800876 SRR31800875 SRR31800874)

download_one() {
    local SRR="$1"
    local LABEL="$2"
    local URL_R1="$3"
    local URL_R2="$4"
    local OUTDIR="$RAW/$LABEL"
    mkdir -p "$OUTDIR"

    for r in 1 2; do
        local FN="${SRR}_${r}.fastq.gz"
        local DEST="$OUTDIR/$FN"
        local URL=$([ "$r" = "1" ] && echo "$URL_R1" || echo "$URL_R2")
        if [ -f "$DEST" ] && [ -s "$DEST" ]; then
            local SIZE=$(stat -c %s "$DEST")
            if [ "$SIZE" -gt 104857600 ]; then
                echo "[$(date -Iseconds)] $LABEL R$r already present (${SIZE} bytes) — skipping."
                continue
            fi
        fi
        echo "[$(date -Iseconds)] $LABEL R$r downloading via lftp pget (8 streams) from $URL ..."
        # lftp pget -n 8 -c uses 8 parallel byte-range streams + resume.
        # Far faster than single-stream wget on HTTP origins that allow ranges (ENA does).
        local URL_DIR=$(dirname "$URL")
        local URL_FN=$(basename "$URL")
        cd "$OUTDIR" && \
        lftp -e "set net:max-retries 20; set net:reconnect-interval-base 30; set xfer:clobber yes; pget -n 8 -c \"$URL\" -o \"$FN\"; bye" \
            || { echo "ERROR: $FN download failed"; return 1; }
        echo "[$(date -Iseconds)] $LABEL R$r DONE ($(stat -c %s "$DEST") bytes)"
    done
}

export -f download_one
export RAW
PIDS=()
for i in "${!SRRS[@]}"; do
    download_one "${SRRS[$i]}" "${LABELS[$i]}" "${URLS_R1[$i]}" "${URLS_R2[$i]}" &
    PIDS+=($!)
done

# Wait for all + capture failures
FAIL=0
for pid in "${PIDS[@]}"; do
    wait "$pid" || FAIL=$((FAIL+1))
done

echo "[$(date -Iseconds)] All download tasks finished. Failures: $FAIL"
ls -lh "$RAW"/*/*.fastq.gz | head -20

# SLURM manifest
PROJ="/project/def-srussell/nirosha_xeno"
MANIFEST="$PROJ/scorial/01_state/slurm_manifests/job_${SLURM_JOB_ID}.yaml"
cat > "$MANIFEST" <<EOF
---
job_id: ${SLURM_JOB_ID}
job_name: hA_walker_dl
task_id: T-09a-dl
stage: comparator_walker_download
purpose: "DEC-017: download raw 10X scRNA FASTQs from ENA for PRJNA1200678 (Walker 2025 axolotl brain regen). 4 samples, ~485 GB compressed FASTQ to /scratch/srussell/srussell/nirosha_xeno_comparators/walker/raw/."
end_time: $(date -Iseconds)
node: $(hostname)
exit_code: $FAIL
samples:
  - SRR31800877 / uninjured_brain  (KE1, 417M reads)
  - SRR31800876 / tail_brain       (KE2, 420M reads)
  - SRR31800875 / limb_brain       (372M reads)
  - SRR31800874 / muscimol_brain   (407M reads)
outputs:
  raw_dir: $RAW
EOF
echo "Manifest: $MANIFEST"
exit $FAIL
