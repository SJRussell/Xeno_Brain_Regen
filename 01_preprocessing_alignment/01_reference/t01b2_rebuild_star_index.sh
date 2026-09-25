#!/bin/bash
#SBATCH --job-name=t01b2_star_idx_v2
#SBATCH --account=def-srussell
#SBATCH --time=01:30:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=16
#SBATCH --output=/scratch/srussell/srussell/logs/t01b2_star_idx_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/t01b2_star_idx_%j.err

set -euo pipefail

module purge
module load StdEnv/2023 star/2.7.11b

PROJECT_ROOT=/project/def-srussell/nirosha_xeno
REF_DIR="${PROJECT_ROOT}/ref/XENLA_10.1"
INDEX_DIR="${REF_DIR}/star_index"
MANIFEST_DIR="${PROJECT_ROOT}/scorial/01_state/slurm_manifests"
mkdir -p "$MANIFEST_DIR"

cd "$REF_DIR"

# REBUILD: previous index built with default --sjdbGTFtagExonParentGene
# (=gene_id, GTF convention) which doesn't exist on NCBI GFF3 exon lines.
# NCBI GFF3 exon lines carry gene=LOC<NNN>; that attribute name works.
# Add --sjdbGTFtagExonParentGeneName gene for symbol-level output too.
SJDB_OVERHANG=99

# Wipe broken index dir
rm -rf "${INDEX_DIR}"
mkdir -p "${INDEX_DIR}"

echo "[$(date -Iseconds)] STAR genomeGenerate (v2 with proper GFF3 gene-id mapping)"

STAR \
    --runMode genomeGenerate \
    --genomeDir "${INDEX_DIR}" \
    --genomeFastaFiles "${REF_DIR}/genome.fa" \
    --sjdbGTFfile "${REF_DIR}/annotation_GCF.gff3" \
    --sjdbGTFtagExonParentTranscript Parent \
    --sjdbGTFtagExonParentGene gene \
    --sjdbGTFtagExonParentGeneName gene \
    --sjdbGTFtagExonParentGeneType gene_biotype \
    --sjdbOverhang ${SJDB_OVERHANG} \
    --genomeSAindexNbases 14 \
    --runThreadN ${SLURM_CPUS_PER_TASK}

EXIT_CODE=$?
echo "[$(date -Iseconds)] STAR exited ${EXIT_CODE}"

ls -la "${INDEX_DIR}"

# Verification: geneInfo.tab should now have ~30-45k entries
GENE_COUNT=$(wc -l < "${INDEX_DIR}/geneInfo.tab")
TX_COUNT=$(wc -l < "${INDEX_DIR}/transcriptInfo.tab")
SJDB_COUNT=$(wc -l < "${INDEX_DIR}/sjdbList.fromGTF.out.tab")
echo "geneInfo.tab line count:           ${GENE_COUNT}"
echo "transcriptInfo.tab line count:      ${TX_COUNT}"
echo "sjdbList.fromGTF.out.tab line count: ${SJDB_COUNT}"

# Sample first 5 gene entries
echo "First 5 entries of geneInfo.tab:"
head -5 "${INDEX_DIR}/geneInfo.tab"

MAN="${MANIFEST_DIR}/job_${SLURM_JOB_ID}.yaml"
{
    echo "---"
    echo "job_id: ${SLURM_JOB_ID}"
    echo "job_name: ${SLURM_JOB_NAME}"
    echo "task_id: T-01b (rebuild)"
    echo "stage: stage_ref_xenla_10.1_star_index_v2"
    echo "exit_code: ${EXIT_CODE}"
    echo "end_time: $(date -Iseconds)"
    echo "node: $(hostname)"
    echo "tool: STAR/2.7.11b"
    echo "supersedes_job: 13047888"
    echo "fix_applied: 'Added --sjdbGTFtagExonParentGene gene (NCBI GFF3 attribute) to populate geneInfo.tab; v1 index produced single MissingGeneID entry.'"
    echo "parameters:"
    echo "  sjdbOverhang: ${SJDB_OVERHANG}"
    echo "  sjdbGTFtagExonParentTranscript: Parent"
    echo "  sjdbGTFtagExonParentGene: gene"
    echo "  sjdbGTFtagExonParentGeneName: gene"
    echo "  sjdbGTFtagExonParentGeneType: gene_biotype"
    echo "  threads: ${SLURM_CPUS_PER_TASK}"
    echo "verification:"
    echo "  geneInfo_lines: ${GENE_COUNT}"
    echo "  transcriptInfo_lines: ${TX_COUNT}"
    echo "  sjdb_from_gtf_lines: ${SJDB_COUNT}"
    echo "output_files:"
    if [[ -d "${INDEX_DIR}" ]]; then
        for f in SA SAindex Genome geneInfo.tab transcriptInfo.tab sjdbList.out.tab chrNameLength.txt exonInfo.tab; do
            if [[ -f "${INDEX_DIR}/${f}" ]]; then
                echo "  - path: ${INDEX_DIR}/${f}"
                echo "    size_bytes: $(stat -c%s "${INDEX_DIR}/${f}")"
            fi
        done
    fi
} > "$MAN"
cat "$MAN"

echo "[$(date -Iseconds)] T-01b (rebuild) complete"
