#!/bin/bash
#SBATCH --job-name=v7_le_extract
#SBATCH --account=def-srussell
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=/scratch/srussell/srussell/logs/v7_le_extract_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/v7_le_extract_%j.err

set -euo pipefail
module purge
module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21
export R_LIBS_USER=/home/srussell/Rlibs/4.5

WK=/scratch/srussell/srussell
PROJ=/home/srussell/projects/def-srussell/nirosha_xeno
SCRIPT=$WK/scripts/extract_leading_edge_v7.R
LOG=$WK/logs
MANIFEST_DIR=$PROJ/scorial/01_state/slurm_manifests
mkdir -p "$LOG" "$MANIFEST_DIR"

cd "$PROJ"
START_ISO=$(date -Iseconds)
echo "[$START_ISO] START v7 leading-edge extraction (job $SLURM_JOB_ID)"

Rscript "$SCRIPT"

END_ISO=$(date -Iseconds)
echo "[$END_ISO] END   v7 leading-edge extraction"

# Manifest
cat > "$MANIFEST_DIR/job_${SLURM_JOB_ID}.yaml" <<EOF
job_id: $SLURM_JOB_ID
job_name: v7_le_extract
purpose: |
  Extract leading-edge gene sets for v7 figure heatmaps from
  pathway/homeolog_sum_strand_corrected_results.rds.
  Targets: Fig 2B (BLAST 7DPA EMT), Fig 3B (BRAIN 24HPA TNFA_NFKB),
           Fig 4B (BRAIN 7DPA ALLOGRAFT_REJECTION),
           Fig 4C (BRAIN 7DPA HEME_METABOLISM).
script_R: $SCRIPT
script_sh: $WK/scripts/extract_leading_edge_v7.sh
input:
  - pathway/homeolog_sum_strand_corrected_results.rds
output:
  - pathway/leading_edge_v7.tsv
walltime_requested: "00:10:00"
cpus: 1
mem: 8G
log_out: $LOG/v7_le_extract_${SLURM_JOB_ID}.out
log_err: $LOG/v7_le_extract_${SLURM_JOB_ID}.err
start: $START_ISO
end: $END_ISO
state: COMPLETED
EOF
echo "Manifest written: $MANIFEST_DIR/job_${SLURM_JOB_ID}.yaml"
