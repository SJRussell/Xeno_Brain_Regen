#!/bin/bash
#SBATCH --job-name=hC_erdei
#SBATCH --account=def-srussell
#SBATCH --time=04:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=8
#SBATCH --output=/scratch/srussell/srussell/logs/hC_erdei_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/hC_erdei_%j.err

# H-C Erdei comparator (DEC-016): mouse cortex bulk RNA-seq from PRJNA1137110.
# 3 cortex samples: SRX25371062 CTX_nt, SRX25371064 CTX_bi, SRX25371066 CTX_fa.
# Pipeline: prefetch + fasterq-dump -> trimmomatic SE -> salmon quant
#   (decoy-free cdna-only) -> tximport (countsFromAbundance=lengthScaledTPM)
#   -> limma-voom -> Hallmark/GO:BP fgsea -> cross-species concordance vs
#   Xenopus 24HPA CT brain (matches Erdei's 1h post-injury timepoint).
# Mirrors hB_descalzi but with quantification (Descalzi had pre-computed RPKM).

set -euo pipefail
module purge
module load StdEnv/2023 sra-toolkit/3.0.9 salmon/1.10.2 r/4.5.0 r-bundle-bioconductor/3.21
export R_LIBS_USER=/home/srussell/Rlibs/4.5

PROJ="/project/def-srussell/nirosha_xeno"
WORKROOT="/scratch/srussell/srussell/nirosha_xeno_comparators/erdei"
RAW="$WORKROOT/raw"
WORK="$WORKROOT/work"
QUANT="$WORK/salmon"
REF_DIR="/scratch/srussell/srussell/refs/mm39_gencode_m33"
mkdir -p "$RAW" "$WORK" "$QUANT" "$REF_DIR"

cd "$WORKROOT"

# ------------------------------------------------------------
# 1. Salmon index for GENCODE mouse M33 (cdna only, decoy-free)
# ------------------------------------------------------------
SALMON_IDX="$REF_DIR/salmon_idx"
TX_FA="$REF_DIR/gencode.vM33.transcripts.fa.gz"
TX2GENE="$REF_DIR/tx2gene.gencode.vM33.tsv"

if [ ! -d "$SALMON_IDX" ] || [ ! -f "$SALMON_IDX/info.json" ]; then
    echo "[$(date -Iseconds)] Building salmon index for GENCODE M33..."
    if [ ! -f "$TX_FA" ]; then
        wget -q -O "$TX_FA" \
          "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/gencode.vM33.transcripts.fa.gz"
    fi
    salmon index -t "$TX_FA" -i "$SALMON_IDX" -k 31 --threads "$SLURM_CPUS_PER_TASK" --keepDuplicates
    echo "[$(date -Iseconds)] Salmon index built."
else
    echo "[$(date -Iseconds)] Salmon index cached at $SALMON_IDX — skipping build."
fi

# Build tx2gene from the FASTA headers (GENCODE format: >ENSMUST...|ENSMUSG...|...|gene_symbol|...|...)
if [ ! -f "$TX2GENE" ]; then
    echo "[$(date -Iseconds)] Building tx2gene from GENCODE FASTA headers..."
    zcat "$TX_FA" | awk -F'|' '/^>/ {sub(/^>/,"",$1); print $1"\t"$2"\t"$6}' > "$TX2GENE"
    echo "tx2gene rows: $(wc -l < "$TX2GENE")"
fi

# ------------------------------------------------------------
# 2. SRA download + fasterq-dump (3 cortex samples)
# ------------------------------------------------------------
# DEC-016 cortex subset SRR mapping (resolved from SRX via ENA filereport):
#   SRR29873993  CTX_nt  (SRX25371062, 17.3M reads)
#   SRR29873991  CTX_bi  (SRX25371064, 16.5M reads)
#   SRR29873989  CTX_fa  (SRX25371066, 20.1M reads)
# fasterq-dump does not accept SRX (container) accessions, only SRR.
SAMPLES=(SRR29873993 SRR29873991 SRR29873989)
NAMES=(CTX_nt CTX_bi CTX_fa)

cd "$RAW"
for i in "${!SAMPLES[@]}"; do
    SRR="${SAMPLES[$i]}"
    NAME="${NAMES[$i]}"
    if [ ! -d "$SRR" ] || [ ! -f "$RAW/$SRR/$SRR.sra" ]; then
        echo "[$(date -Iseconds)] prefetch $SRR..."
        prefetch --max-size 50g "$SRR"
    else
        echo "[$(date -Iseconds)] SRA cached for $SRR — skipping prefetch."
    fi
    if [ ! -f "$RAW/${NAME}.fastq.gz" ]; then
        echo "[$(date -Iseconds)] fasterq-dump $SRR (single-end mouse RNA-seq)..."
        # Erdei is SE Illumina — fasterq-dump default behavior dumps ${SRR}.fastq for SE
        fasterq-dump --threads "$SLURM_CPUS_PER_TASK" --skip-technical \
            --outdir "$RAW" "$RAW/$SRR/$SRR.sra"
        if [ -f "$RAW/${SRR}.fastq" ]; then
            mv "$RAW/${SRR}.fastq" "$RAW/${NAME}.fastq" && gzip "$RAW/${NAME}.fastq"
        elif [ -f "$RAW/${SRR}_1.fastq" ] && [ ! -s "$RAW/${SRR}_2.fastq" ]; then
            mv "$RAW/${SRR}_1.fastq" "$RAW/${NAME}.fastq" && gzip "$RAW/${NAME}.fastq"
            rm -f "$RAW/${SRR}_2.fastq"
        elif [ -f "$RAW/${SRR}_1.fastq" ] && [ -s "$RAW/${SRR}_2.fastq" ]; then
            # Unexpected paired output — concat as single-end (paired R1+R2 treated as more reads)
            cat "$RAW/${SRR}_1.fastq" "$RAW/${SRR}_2.fastq" > "$RAW/${NAME}.fastq" && gzip "$RAW/${NAME}.fastq"
            rm -f "$RAW/${SRR}_1.fastq" "$RAW/${SRR}_2.fastq"
        else
            echo "ERROR: fasterq-dump produced no FASTQ for $SRR"; exit 64
        fi
    else
        echo "[$(date -Iseconds)] FASTQ for $NAME already present — skipping."
    fi
done

ls -lh "$RAW/"*.fastq.gz || true

# ------------------------------------------------------------
# 3. Salmon quant per sample
# ------------------------------------------------------------
for NAME in "${NAMES[@]}"; do
    OUT="$QUANT/$NAME"
    if [ ! -f "$OUT/quant.sf" ]; then
        echo "[$(date -Iseconds)] salmon quant $NAME..."
        FQ="$RAW/${NAME}.fastq.gz"
        if [ ! -f "$FQ" ]; then
            echo "ERROR: missing FASTQ $FQ — aborting"
            exit 1
        fi
        salmon quant -i "$SALMON_IDX" -l A \
            -r "$FQ" \
            -p "$SLURM_CPUS_PER_TASK" \
            --validateMappings --gcBias --seqBias \
            -o "$OUT"
    else
        echo "salmon quant for $NAME already done — skipping."
    fi
done

# ------------------------------------------------------------
# 4. R analysis (tximport + voom + GSEA + concordance)
# ------------------------------------------------------------
echo "[$(date -Iseconds)] Running R analysis..."
Rscript /scratch/srussell/srussell/scripts/hC_erdei.R

# ------------------------------------------------------------
# 5. SLURM manifest (per CR-* manifest spec)
# ------------------------------------------------------------
MANIFEST="$PROJ/scorial/01_state/slurm_manifests/job_${SLURM_JOB_ID}.yaml"
cat > "$MANIFEST" <<EOF
---
job_id: ${SLURM_JOB_ID}
job_name: hC_erdei
task_id: T-09c
stage: comparator_erdei
purpose: "Erdei 2025 IJMS PRJNA1137110 mouse S1 cortex 1h-post-injury bulk RNA-seq comparator. salmon -> tximport (lengthScaledTPM) -> voom -> limma -> Hallmark/GOBP fgsea. Cross-species concordance vs Xenopus 24HPA CT brain. DEC-016 (re-engaged from erroneous drop)."
end_time: $(date -Iseconds)
node: $(hostname)
exit_code: 0
samples:
  - SRX25371062 / CTX_nt (non-treated cortex)
  - SRX25371064 / CTX_bi (burn-injury cortex)
  - SRX25371066 / CTX_fa (formalin cortex)
contrast_primary: bi_vs_nt
sensitivity: include_fa_as_secondary_contrast
outputs:
  workspace: $WORKROOT
  de: /project/def-srussell/nirosha_xeno/comparators/erdei/de_*.tsv
  gsea: /project/def-srussell/nirosha_xeno/comparators/erdei/gsea_*.tsv
  concordance: /project/def-srussell/nirosha_xeno/comparators/erdei/xenopus_concordance.tsv
EOF
echo "Manifest: $MANIFEST"
echo "[$(date -Iseconds)] hC_erdei.sh complete."
