#!/bin/bash
#SBATCH --job-name=hD_liao_dl
#SBATCH --account=def-srussell
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --output=/scratch/srussell/srussell/logs/hD_liao_dl_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/hD_liao_dl_%j.err

set -euo pipefail

RAW=/scratch/srussell/srussell/nirosha_xeno_comparators/liao/raw
mkdir -p "$RAW"
cd "$RAW"

# dge_cell_info.zip (16.4 MB) — full XCL atlas cell metadata
[ -f dge_cell_info.zip ] || curl -sL --max-time 600 -o dge_cell_info.zip \
    "https://ndownloader.figshare.com/files/34026644"

# dge_raw_data.zip (909.7 MB) — raw count matrices per sample
[ -f dge_raw_data.zip ] || curl -sL --max-time 1800 -o dge_raw_data.zip \
    "https://ndownloader.figshare.com/files/34026647"

ls -lh dge_*.zip
echo "---"
# Inspect zip contents (file listing only — no extraction yet)
echo "=== dge_cell_info.zip listing ==="
unzip -l dge_cell_info.zip | head -20
echo "..."
echo "=== dge_raw_data.zip listing (first 30 entries + count) ==="
unzip -l dge_raw_data.zip | head -30
echo "Total entries:"
unzip -l dge_raw_data.zip | tail -3

# Extract dge_cell_info.zip immediately (small)
unzip -o dge_cell_info.zip -d cell_info/
ls cell_info/ | head

# For raw data, only extract Brain* files now (the ZIP probably has per-sample files;
# we only need the Brain ones for our analysis to save disk + time)
echo "=== Brain* entries in raw zip ==="
unzip -l dge_raw_data.zip | grep -i brain | head -10

# Manifest
MANIFEST="/project/def-srussell/nirosha_xeno/scorial/01_state/slurm_manifests/job_${SLURM_JOB_ID}.yaml"
cat > "$MANIFEST" <<EOF
---
job_id: ${SLURM_JOB_ID}
job_name: hD_liao_dl
task_id: T-09d-download
stage: comparator_liao_setup
purpose: "Download Liao 2022 XCL processed-data Figshare deposit (cell_info.zip + dge_raw_data.zip). Extract cell_info; survey raw data zip for brain-tissue files."
end_time: $(date -Iseconds)
node: $(hostname)
exit_code: 0
outputs:
  raw_dir: /scratch/srussell/srussell/nirosha_xeno_comparators/liao/raw
  cell_info_csv: /scratch/srussell/srussell/nirosha_xeno_comparators/liao/raw/Xenopus_Figure1_cell_info.csv
  cell_info_zip_extracted: /scratch/srussell/srussell/nirosha_xeno_comparators/liao/raw/cell_info/
  raw_data_zip_listing: see job stdout
EOF
echo "Manifest written."
