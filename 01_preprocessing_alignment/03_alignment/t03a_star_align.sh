#!/bin/bash
#SBATCH --job-name=t03a_star_align
#SBATCH --account=def-srussell
#SBATCH --time=00:45:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-30%8
#SBATCH --output=/scratch/srussell/srussell/logs/t03a_star_%A_%a.out
#SBATCH --error=/scratch/srussell/srussell/logs/t03a_star_%A_%a.err

set -euo pipefail

module purge
module load StdEnv/2023 star/2.7.11b samtools/1.20

PROJECT_ROOT=/project/def-srussell/nirosha_xeno
TRIM_DIR="${PROJECT_ROOT}/trim"
STAR_INDEX="${PROJECT_ROOT}/ref/XENLA_10.1/star_index"
STAR_OUT="${PROJECT_ROOT}/star"
MANIFEST_DIR="${PROJECT_ROOT}/scorial/01_state/slurm_manifests"
mkdir -p "$STAR_OUT" "$MANIFEST_DIR"

# Pick the SLURM_ARRAY_TASK_ID-th trimmed FASTQ
FASTQS=( "${TRIM_DIR}"/*.trim.fastq.gz )
IDX=$((SLURM_ARRAY_TASK_ID - 1))
if [[ $IDX -ge ${#FASTQS[@]} ]]; then
    echo "ERROR: array index $IDX out of range (have ${#FASTQS[@]} samples)"
    exit 2
fi

IN_FQ="${FASTQS[$IDX]}"
SAMPLE=$(basename "$IN_FQ" .trim.fastq.gz)
SAMPLE_OUT="${STAR_OUT}/${SAMPLE}"
mkdir -p "$SAMPLE_OUT"

echo "[$(date -Iseconds)] STAR align: $SAMPLE"
echo "  input: $IN_FQ"
echo "  output: $SAMPLE_OUT"

cd "$SAMPLE_OUT"

# Idempotent skip if BAM exists and is non-empty
if [[ -s "${SAMPLE_OUT}/Aligned.sortedByCoord.out.bam" ]]; then
    echo "  [SKIP] BAM already exists and is non-empty"
else
    STAR \
        --runMode alignReads \
        --genomeDir "$STAR_INDEX" \
        --readFilesIn "$IN_FQ" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outFileNamePrefix "${SAMPLE_OUT}/" \
        --quantMode GeneCounts \
        --twopassMode Basic \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --runThreadN ${SLURM_CPUS_PER_TASK}
fi

# Index BAM for downstream tools
samtools index -@ ${SLURM_CPUS_PER_TASK} "${SAMPLE_OUT}/Aligned.sortedByCoord.out.bam"

EXIT_CODE=$?

# Report
N_INPUT=$(grep "Number of input reads" "${SAMPLE_OUT}/Log.final.out" | awk '{print $NF}')
PCT_UNIQUE=$(grep "Uniquely mapped reads %" "${SAMPLE_OUT}/Log.final.out" | awk '{print $NF}')
PCT_MULTI=$(grep "% of reads mapped to multiple loci" "${SAMPLE_OUT}/Log.final.out" | awk '{print $NF}')

echo "[$(date -Iseconds)] $SAMPLE done: input=${N_INPUT}, unique=${PCT_UNIQUE}, multi=${PCT_MULTI}"

MAN="${MANIFEST_DIR}/job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.yaml"
{
    echo "---"
    echo "job_id: ${SLURM_JOB_ID}"
    echo "array_id: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
    echo "job_name: ${SLURM_JOB_NAME}"
    echo "task_id: T-03a"
    echo "stage: stage_B_star_align"
    echo "sample: ${SAMPLE}"
    echo "exit_code: ${EXIT_CODE}"
    echo "end_time: $(date -Iseconds)"
    echo "node: $(hostname)"
    echo "tool: STAR/2.7.11b"
    echo "metrics:"
    echo "  input_reads: ${N_INPUT}"
    echo "  uniquely_mapped_pct: ${PCT_UNIQUE}"
    echo "  multimapped_pct: ${PCT_MULTI}"
    echo "outputs:"
    echo "  bam: ${SAMPLE_OUT}/Aligned.sortedByCoord.out.bam"
    echo "  log: ${SAMPLE_OUT}/Log.final.out"
    echo "  reads_per_gene: ${SAMPLE_OUT}/ReadsPerGene.out.tab"
    echo "  splice_junctions: ${SAMPLE_OUT}/SJ.out.tab"
} > "$MAN"

echo "[$(date -Iseconds)] T-03a sample $SAMPLE complete"
