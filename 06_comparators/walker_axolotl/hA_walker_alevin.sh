#!/bin/bash
#SBATCH --job-name=hA_walker_alevin
#SBATCH --account=def-srussell
#SBATCH --time=04:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-3
#SBATCH --output=/scratch/srussell/srussell/logs/hA_walker_alevin_%A_%a.out
#SBATCH --error=/scratch/srussell/srussell/logs/hA_walker_alevin_%A_%a.err

# H-A Walker alevin quantification (DEC-017): salmon alevin per-sample
# 10X Chromium v3 quant. Submit with --dependency=afterok:<DL>:<REF> so it
# auto-launches when both the FASTQ download (13126403) and reference build
# (13129967) complete.
#
# Output: per-sample alevin quants_mat.gz + permit/raw barcode freq
#         under /scratch/srussell/srussell/nirosha_xeno_comparators/walker/work/alevin/<LABEL>/

set -euo pipefail
module purge
module load StdEnv/2023 salmon/1.10.2

WORKROOT="/scratch/srussell/srussell/nirosha_xeno_comparators/walker"
RAW="$WORKROOT/raw"
WORK="$WORKROOT/work"
ALEVIN="$WORK/alevin"
mkdir -p "$ALEVIN"

REF_DIR="/scratch/srussell/srussell/refs/amext_v47"
SALMON_IDX="$REF_DIR/salmon_idx"
TX2GENE="$REF_DIR/tx2gene_2col.tsv"  # 2-column form (tx_id, gene_id) — alevin standard

# Per-task sample assignment (array index 0..3)
LABELS=(uninjured_brain tail_brain limb_brain muscimol_brain)
SRRS=(SRR31800877 SRR31800876 SRR31800875 SRR31800874)
LABEL="${LABELS[$SLURM_ARRAY_TASK_ID]}"
SRR="${SRRS[$SLURM_ARRAY_TASK_ID]}"

R1="$RAW/$LABEL/${SRR}_1.fastq.gz"
R2="$RAW/$LABEL/${SRR}_2.fastq.gz"
OUTDIR="$ALEVIN/$LABEL"

echo "[$(date -Iseconds)] alevin quant: $LABEL ($SRR)"
echo "  R1: $R1 ($(stat -c %s "$R1" 2>/dev/null || echo MISSING) bytes)"
echo "  R2: $R2 ($(stat -c %s "$R2" 2>/dev/null || echo MISSING) bytes)"

[ -f "$R1" ] || { echo "ERROR: missing R1"; exit 4; }
[ -f "$R2" ] || { echo "ERROR: missing R2"; exit 4; }
[ -f "$SALMON_IDX/info.json" ] || { echo "ERROR: salmon index not ready"; exit 5; }

# salmon alevin handles 10X Chromium v3 (CB=16, UMI=12) via --chromiumV3
# --tgMap maps transcripts -> genes (tx2gene first 2 cols)
salmon alevin -l ISR \
    -1 "$R1" \
    -2 "$R2" \
    -i "$SALMON_IDX" \
    --chromiumV3 \
    -p "$SLURM_CPUS_PER_TASK" \
    --tgMap "$TX2GENE" \
    --dumpFeatures \
    -o "$OUTDIR"

echo "[$(date -Iseconds)] alevin done for $LABEL"
ls -la "$OUTDIR" | head -10

# SLURM manifest (per-task)
PROJ="/project/def-srussell/nirosha_xeno"
MANIFEST="$PROJ/scorial/01_state/slurm_manifests/job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.yaml"
cat > "$MANIFEST" <<EOF
---
job_id: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
job_name: hA_walker_alevin
task_id: T-09a-quant
stage: comparator_walker_alevin
sample_label: $LABEL
sra_run: $SRR
purpose: "salmon alevin quant of Walker 10X scRNA sample $LABEL against AmexT_v47 transcriptome."
end_time: $(date -Iseconds)
node: $(hostname)
exit_code: 0
outputs:
  - $OUTDIR/alevin/quants_mat.gz
EOF
echo "Manifest: $MANIFEST"
