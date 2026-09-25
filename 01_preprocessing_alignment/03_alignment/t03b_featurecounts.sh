#!/bin/bash
#SBATCH --job-name=t03b_fc
#SBATCH --account=def-srussell
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --output=/scratch/srussell/srussell/logs/t03b_fc_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/t03b_fc_%j.err

set -euo pipefail

module purge
module load StdEnv/2023 subread/2.0.6

PROJECT_ROOT=/project/def-srussell/nirosha_xeno
STAR_OUT="${PROJECT_ROOT}/star"
COUNTS_DIR="${PROJECT_ROOT}/counts"
GFF3="${PROJECT_ROOT}/ref/XENLA_10.1/annotation_GCF.gff3"
MANIFEST_DIR="${PROJECT_ROOT}/scorial/01_state/slurm_manifests"
mkdir -p "$COUNTS_DIR" "$MANIFEST_DIR"

# Collect BAMs in stable sample-name order
BAMS=( "${STAR_OUT}"/*/Aligned.sortedByCoord.out.bam )
echo "[$(date -Iseconds)] featureCounts on ${#BAMS[@]} BAMs"
echo "  GFF3: $GFF3"

# featureCounts on GFF3 — use -F GFF, -t exon, -g gene (NCBI GFF3 attribute), -O for overlap.
# -s 2 = reverse-stranded (matches BBcore 2019 preprocess; Salmon also uses SR libtype).
# -p NOT used (single-end reads).

featureCounts \
    -T ${SLURM_CPUS_PER_TASK} \
    -F GFF \
    -t exon \
    -g gene \
    -O \
    -s 2 \
    -a "$GFF3" \
    -o "${COUNTS_DIR}/fc_counts.gene.txt" \
    "${BAMS[@]}"

EXIT_CODE=$?

# Verify
N_GENES=$(awk 'NR>2 {n++} END{print n}' "${COUNTS_DIR}/fc_counts.gene.txt")
N_SAMPLES_OUT=$(head -2 "${COUNTS_DIR}/fc_counts.gene.txt" | tail -1 | awk '{print NF-6}')
echo "[$(date -Iseconds)] featureCounts done: ${N_GENES} genes × ${N_SAMPLES_OUT} samples"

# Sanity — assignment percentages from .summary file
echo "=== featureCounts summary (assignment rates) ==="
cat "${COUNTS_DIR}/fc_counts.gene.txt.summary"

MAN="${MANIFEST_DIR}/job_${SLURM_JOB_ID}.yaml"
{
    echo "---"
    echo "job_id: ${SLURM_JOB_ID}"
    echo "job_name: ${SLURM_JOB_NAME}"
    echo "task_id: T-03b"
    echo "stage: stage_B_featurecounts"
    echo "exit_code: ${EXIT_CODE}"
    echo "end_time: $(date -Iseconds)"
    echo "node: $(hostname)"
    echo "tool: subread/2.0.6 featureCounts"
    echo "parameters:"
    echo "  -t exon -g gene -O -s 2 -F GFF"
    echo "n_genes: ${N_GENES}"
    echo "n_samples: ${N_SAMPLES_OUT}"
    echo "outputs:"
    echo "  counts_table: ${COUNTS_DIR}/fc_counts.gene.txt"
    echo "  summary:      ${COUNTS_DIR}/fc_counts.gene.txt.summary"
} > "$MAN"

echo "[$(date -Iseconds)] T-03b complete"
