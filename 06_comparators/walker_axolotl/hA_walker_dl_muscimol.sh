#!/bin/bash
#SBATCH --job-name=walker_dl_musc
#SBATCH --account=def-srussell
#SBATCH --time=06:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4
#SBATCH --output=/scratch/srussell/srussell/logs/walker_dl_musc_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/walker_dl_musc_%j.err

set -euo pipefail
RAW="/scratch/srussell/srussell/nirosha_xeno_comparators/walker/raw/muscimol_brain"
URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR318/074/SRR31800874/SRR31800874_2.fastq.gz"
DEST="$RAW/SRR31800874_2.fastq.gz"
cd "$RAW"
echo "[$(date -Iseconds)] Re-downloading muscimol R2..."
lftp -e "set net:max-retries 30; pget -n 8 -c \"$URL\" -o $(basename "$DEST"); bye"
echo "[$(date -Iseconds)] Validating gzip..."
gzip -t "$DEST" && echo "OK" || { echo "STILL CORRUPT"; exit 1; }
