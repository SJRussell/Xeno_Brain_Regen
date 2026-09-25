#!/bin/bash
#SBATCH --job-name=t03c_salmon_quant
#SBATCH --account=def-srussell
#SBATCH --time=00:30:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-30%8
#SBATCH --output=/scratch/srussell/srussell/logs/t03c_salmon_%A_%a.out
#SBATCH --error=/scratch/srussell/srussell/logs/t03c_salmon_%A_%a.err

set -euo pipefail

module purge
module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 salmon/1.10.2

PROJECT_ROOT=/project/def-srussell/nirosha_xeno
TRIM_DIR="${PROJECT_ROOT}/trim"
SALMON_INDEX="${PROJECT_ROOT}/ref/XENLA_10.1/salmon_index"
SALMON_OUT="${PROJECT_ROOT}/salmon"
MANIFEST_DIR="${PROJECT_ROOT}/scorial/01_state/slurm_manifests"
mkdir -p "$SALMON_OUT" "$MANIFEST_DIR"

FASTQS=( "${TRIM_DIR}"/*.trim.fastq.gz )
IDX=$((SLURM_ARRAY_TASK_ID - 1))
if [[ $IDX -ge ${#FASTQS[@]} ]]; then
    echo "ERROR: array index $IDX out of range"
    exit 2
fi

IN_FQ="${FASTQS[$IDX]}"
SAMPLE=$(basename "$IN_FQ" .trim.fastq.gz)
SAMPLE_OUT="${SALMON_OUT}/${SAMPLE}"

echo "[$(date -Iseconds)] Salmon quant: $SAMPLE"

# Skip if already done
if [[ -s "${SAMPLE_OUT}/quant.sf" ]]; then
    echo "  [SKIP] quant.sf already exists"
else
    # libtype SR = stranded reverse, single-end (matches featureCounts -s 2)
    salmon quant \
        -i "$SALMON_INDEX" \
        -l SR \
        -r "$IN_FQ" \
        --validateMappings \
        --gcBias \
        --seqBias \
        -p ${SLURM_CPUS_PER_TASK} \
        -o "$SAMPLE_OUT"
fi

EXIT_CODE=$?

# Read mapping rate from meta_info.json
MAP_RATE=$(python3 -c "
import json
with open('${SAMPLE_OUT}/aux_info/meta_info.json') as f:
    d = json.load(f)
print(f\"{d['percent_mapped']:.2f}\")
" 2>/dev/null || echo "NA")
echo "[$(date -Iseconds)] $SAMPLE done: mapping rate=${MAP_RATE}%"

MAN="${MANIFEST_DIR}/job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_salmon.yaml"
{
    echo "---"
    echo "job_id: ${SLURM_JOB_ID}"
    echo "array_id: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
    echo "task_id: T-03c"
    echo "stage: stage_B_salmon_quant"
    echo "sample: ${SAMPLE}"
    echo "exit_code: ${EXIT_CODE}"
    echo "end_time: $(date -Iseconds)"
    echo "node: $(hostname)"
    echo "tool: salmon/1.10.2"
    echo "parameters:"
    echo "  libtype: SR"
    echo "  validateMappings: true"
    echo "  gcBias: true"
    echo "  seqBias: true"
    echo "metrics:"
    echo "  percent_mapped: ${MAP_RATE}"
    echo "outputs:"
    echo "  quant_sf: ${SAMPLE_OUT}/quant.sf"
} > "$MAN"

echo "[$(date -Iseconds)] T-03c sample $SAMPLE complete"
