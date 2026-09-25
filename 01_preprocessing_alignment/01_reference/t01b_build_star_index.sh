#!/bin/bash
#SBATCH --job-name=t01b_star_idx
#SBATCH --account=def-srussell
#SBATCH --time=06:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=16
#SBATCH --output=/scratch/srussell/srussell/logs/t01b_star_idx_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/t01b_star_idx_%j.err

set -euo pipefail

module purge
module load StdEnv/2023 star/2.7.11b

PROJECT_ROOT=/project/def-srussell/nirosha_xeno
REF_DIR="${PROJECT_ROOT}/ref/XENLA_10.1"
INDEX_DIR="${REF_DIR}/star_index"
MANIFEST_DIR="${PROJECT_ROOT}/scorial/01_state/slurm_manifests"
mkdir -p "$INDEX_DIR" "$MANIFEST_DIR"

cd "$REF_DIR"

# sjdbOverhang chosen as 99 for compatibility with 50-100bp reads.
# Pipeline plan suggests 70; using 99 as a safe permissive default
# until read length is verified from raw FASTQs in T-02.
SJDB_OVERHANG=99

echo "[$(date -Iseconds)] Starting STAR genomeGenerate"
echo "  Genome: ${REF_DIR}/genome.fa"
echo "  GFF3:   ${REF_DIR}/annotation_GCF.gff3"
echo "  sjdbOverhang: ${SJDB_OVERHANG}"

STAR \
    --runMode genomeGenerate \
    --genomeDir "${INDEX_DIR}" \
    --genomeFastaFiles "${REF_DIR}/genome.fa" \
    --sjdbGTFfile "${REF_DIR}/annotation_GCF.gff3" \
    --sjdbGTFtagExonParentTranscript Parent \
    --sjdbOverhang ${SJDB_OVERHANG} \
    --genomeSAindexNbases 14 \
    --runThreadN ${SLURM_CPUS_PER_TASK}

EXIT_CODE=$?
echo "[$(date -Iseconds)] STAR exited ${EXIT_CODE}"

ls -la "${INDEX_DIR}"

MAN="${MANIFEST_DIR}/job_${SLURM_JOB_ID}.yaml"
{
    echo "---"
    echo "job_id: ${SLURM_JOB_ID}"
    echo "job_name: ${SLURM_JOB_NAME}"
    echo "task_id: T-01b"
    echo "stage: stage_ref_xenla_10.1_star_index"
    echo "exit_code: ${EXIT_CODE}"
    echo "end_time: $(date -Iseconds)"
    echo "node: $(hostname)"
    echo "tool: STAR/2.7.11b"
    echo "parameters:"
    echo "  sjdbOverhang: ${SJDB_OVERHANG}"
    echo "  sjdbGTFtagExonParentTranscript: Parent"
    echo "  threads: ${SLURM_CPUS_PER_TASK}"
    echo "output_files:"
    if [[ -d "${INDEX_DIR}" ]]; then
        for f in SA SAindex Genome geneInfo.tab sjdbList.out.tab chrNameLength.txt; do
            if [[ -f "${INDEX_DIR}/${f}" ]]; then
                echo "  - path: ${INDEX_DIR}/${f}"
                echo "    size_bytes: $(stat -c%s "${INDEX_DIR}/${f}")"
            fi
        done
    fi
} > "$MAN"
cat "$MAN"

echo "[$(date -Iseconds)] T-01b complete"
