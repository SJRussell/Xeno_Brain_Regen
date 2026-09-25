#!/bin/bash
#SBATCH --job-name=hA_walker_ref
#SBATCH --account=def-srussell
#SBATCH --time=03:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --output=/scratch/srussell/srussell/logs/hA_walker_ref_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/hA_walker_ref_%j.err

# H-A Walker reference build (DEC-017): axolotl AmexT_v47 transcriptome
# (matches Walker 2025 reference) + salmon decoy-free index for alevin
# 10X Chromium v3 quantification.
# Lighter footprint than STAR genome (8 Gb genome → 64 GB index).
# salmon alevin handles 10X CB+UMI parsing natively via --chromiumV3.
#
# Outputs:
#   /scratch/srussell/srussell/refs/amext_v47/AmexT_v47.transcripts.fa.gz
#   /scratch/srussell/srussell/refs/amext_v47/AmexT_v47-AmexG_v6.0-DD.gtf.gz
#   /scratch/srussell/srussell/refs/amext_v47/tx2gene.tsv
#   /scratch/srussell/srussell/refs/amext_v47/salmon_idx/

set -euo pipefail
module purge
module load StdEnv/2023 salmon/1.10.2

REF_DIR="/scratch/srussell/srussell/refs/amext_v47"
mkdir -p "$REF_DIR"
cd "$REF_DIR"

TX_FA="$REF_DIR/AmexT_v47.transcripts.fa.gz"
GTF="$REF_DIR/AmexT_v47-AmexG_v6.0-DD.gtf.gz"
TX2GENE="$REF_DIR/tx2gene.tsv"
SALMON_IDX="$REF_DIR/salmon_idx"

# ------------------------------------------------------------
# 1. Download GTF (annotation)
# ------------------------------------------------------------
if [ ! -f "$GTF" ]; then
    echo "[$(date -Iseconds)] Downloading AmexT_v47 GTF..."
    wget --tries=10 --continue --timeout=60 --waitretry=30 --no-verbose \
         -O "$GTF" "http://www.axolotl-omics.org/dl/AmexT_v47-AmexG_v6.0-DD.gtf.gz"
fi
echo "GTF size: $(stat -c %s "$GTF") bytes"

# ------------------------------------------------------------
# 2. Download transcriptome FASTA (DNA sequences via API)
# axolotl-omics API returns large response; allow long timeouts.
# ------------------------------------------------------------
if [ ! -f "$TX_FA" ]; then
    echo "[$(date -Iseconds)] Downloading AmexT_v47 transcriptome (API DNA)..."
    # Two URL candidates: API endpoint (preferred) + direct dl/ endpoint (fallback)
    URL_PRIMARY="http://www.axolotl-omics.org/api?method=Assembly.getSequences&assembly=47&type=dna"
    URL_FALLBACK="http://www.axolotl-omics.org/dl/AmexT_v47-AmexG_v6.0-DD.fa.gz"
    wget --tries=10 --continue --timeout=120 --waitretry=60 --no-verbose \
         -O "$TX_FA" "$URL_PRIMARY" \
        || { echo "Primary URL failed; trying fallback dl/..."; \
             wget --tries=10 --continue --timeout=120 --waitretry=60 --no-verbose \
                  -O "$TX_FA" "$URL_FALLBACK"; }
fi
SIZE=$(stat -c %s "$TX_FA")
echo "Transcriptome size: $SIZE bytes"
if [ "$SIZE" -lt 1000000 ]; then
    echo "ERROR: transcriptome too small ($SIZE bytes) — likely an HTML error page"
    head -5 <(zcat "$TX_FA" 2>/dev/null || cat "$TX_FA")
    exit 4
fi

# Validate FASTA: should start with '>'. Using read-from-process-substitution
# avoids pipefail SIGPIPE (caused exit 13 in 13129967; same class as ERR-2026-04-30 zcat-head bug).
read -r FIRST < <(zcat "$TX_FA" 2>/dev/null) || true
if [[ "$FIRST" != ">"* ]]; then
    echo "ERROR: transcriptome doesn't look like FASTA. First line:"
    echo "$FIRST"
    exit 5
fi
echo "FASTA first record header: $FIRST"

# ------------------------------------------------------------
# 3. Build tx2gene from GTF (transcript_id → gene_id → gene_name where available)
# ------------------------------------------------------------
if [ ! -f "$TX2GENE" ]; then
    echo "[$(date -Iseconds)] Building tx2gene from GTF..."
    zcat "$GTF" | awk -F'\t' '$3 == "transcript" {
        match($9, /transcript_id "[^"]*"/);  tx_id = substr($9, RSTART+15, RLENGTH-16);
        match($9, /gene_id "[^"]*"/);        gn_id = substr($9, RSTART+9,  RLENGTH-10);
        match($9, /gene_name "[^"]*"/);
        if (RSTART > 0) gn_nm = substr($9, RSTART+11, RLENGTH-12); else gn_nm = gn_id;
        print tx_id "\t" gn_id "\t" gn_nm;
    }' > "$TX2GENE"
    echo "tx2gene rows: $(wc -l < "$TX2GENE")"
fi
head -3 "$TX2GENE"

# ------------------------------------------------------------
# 4. Build salmon index (decoy-free, k=31)
# ------------------------------------------------------------
if [ ! -d "$SALMON_IDX" ] || [ ! -f "$SALMON_IDX/info.json" ]; then
    echo "[$(date -Iseconds)] Building salmon decoy-free index..."
    salmon index -t "$TX_FA" -i "$SALMON_IDX" -k 31 \
        --threads "$SLURM_CPUS_PER_TASK" --keepDuplicates
    echo "[$(date -Iseconds)] Salmon index built."
else
    echo "Salmon index already cached at $SALMON_IDX"
fi
du -sh "$SALMON_IDX"

# ------------------------------------------------------------
# SLURM manifest
# ------------------------------------------------------------
PROJ="/project/def-srussell/nirosha_xeno"
MANIFEST="$PROJ/scorial/01_state/slurm_manifests/job_${SLURM_JOB_ID}.yaml"
cat > "$MANIFEST" <<EOF
---
job_id: ${SLURM_JOB_ID}
job_name: hA_walker_ref
task_id: T-09a-ref
stage: comparator_walker_ref_build
purpose: "DEC-017: axolotl AmexT_v47 transcriptome + GTF + salmon decoy-free index for 10X Chromium v3 alevin quantification. Reference matches Walker 2025 published reference."
end_time: $(date -Iseconds)
node: $(hostname)
exit_code: 0
outputs:
  - $TX_FA
  - $GTF
  - $TX2GENE
  - $SALMON_IDX
EOF
echo "Manifest: $MANIFEST"
echo "[$(date -Iseconds)] hA_walker_refbuild.sh complete."
