#!/bin/bash
#SBATCH --job-name=walker_alevin_musc
#SBATCH --account=def-srussell
#SBATCH --time=03:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=8
#SBATCH --output=/scratch/srussell/srussell/logs/walker_alevin_musc_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/walker_alevin_musc_%j.err

# Re-run salmon alevin for the muscimol sample only (post-RED-1 recovery).
# Submitted with --dependency=afterok:<RED-1 job> so it auto-launches once
# the recovered FASTQs are staged.
set -euo pipefail
module purge
module load StdEnv/2023 salmon/1.10.2

WORKROOT="/scratch/srussell/srussell/nirosha_xeno_comparators/walker"
RAW="$WORKROOT/raw"
WORK="$WORKROOT/work"
ALEVIN="$WORK/alevin"

REF_DIR="/scratch/srussell/srussell/refs/amext_v47"
SALMON_IDX="$REF_DIR/salmon_idx"
TX2GENE="$REF_DIR/tx2gene_2col.tsv"

LABEL="muscimol_brain"
SRR="SRR31800874"
R1="$RAW/$LABEL/${SRR}_1.fastq.gz"
R2="$RAW/$LABEL/${SRR}_2.fastq.gz"
OUTDIR="$ALEVIN/$LABEL"

# Move any existing failed-output dir aside (preserve evidence)
if [ -d "$OUTDIR" ]; then
    BAK="${OUTDIR}_failed_pre_RED1_$(date +%Y%m%dT%H%M%S)"
    mv "$OUTDIR" "$BAK"
    echo "[$(date -Iseconds)] preserved prior failed output as $BAK"
fi

echo "[$(date -Iseconds)] alevin quant: $LABEL ($SRR) — RED-1 recovery re-run"
echo "  R1: $R1 ($(stat -c %s "$R1" 2>/dev/null || echo MISSING) bytes)"
echo "  R2: $R2 ($(stat -c %s "$R2" 2>/dev/null || echo MISSING) bytes)"

[ -f "$R1" ] || { echo "ERROR: missing R1 (RED-1 recovery did not produce it)"; exit 4; }
[ -f "$R2" ] || { echo "ERROR: missing R2 (RED-1 recovery did not produce it)"; exit 4; }
[ -f "$SALMON_IDX/info.json" ] || { echo "ERROR: salmon index not ready"; exit 5; }

# Validate gzip integrity before alevin
gzip -t "$R1" || { echo "ERROR: R1 corrupt"; exit 6; }
gzip -t "$R2" || { echo "ERROR: R2 corrupt"; exit 6; }

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
ls -la "$OUTDIR/alevin" | head -10

PROJ="/project/def-srussell/nirosha_xeno"
MANIFEST="$PROJ/scorial/01_state/slurm_manifests/job_${SLURM_JOB_ID}.yaml"
cat > "$MANIFEST" <<EOF
---
job_id: ${SLURM_JOB_ID}
job_name: walker_alevin_musc
task_id: RED-1-alevin
stage: comparator_walker_alevin_muscimol_recovery
purpose: "RED-1 follow-up: re-run salmon alevin on recovered SRR31800874 muscimol FASTQs (post-NCBI-prefetch). Reverses DEC-019 by integrating the muscimol arm."
end_time: $(date -Iseconds)
node: $(hostname)
exit_code: 0
related_decisions: [DEC-017, DEC-019]
inputs:
  - $R1
  - $R2
outputs:
  - $OUTDIR/alevin/quants_mat.gz
EOF
echo "Manifest: $MANIFEST"
