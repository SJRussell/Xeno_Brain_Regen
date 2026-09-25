#!/bin/bash
#SBATCH --job-name=setup_install_libs
#SBATCH --account=def-srussell
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=/scratch/srussell/srussell/logs/setup_install_libs_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/setup_install_libs_%j.err

set -euo pipefail
module purge
module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21
export R_LIBS_USER=/home/srussell/Rlibs/4.5
mkdir -p "$R_LIBS_USER"

# Install MOFA2 + UCell + BayesPrism into user lib.
# These are the only packages missing from the bundle for Stage G + comparators.
# mofapy2 (Python backend) is installed via MOFA2::install_mofapy2() into a per-user
# virtualenv after R-level install.
Rscript - <<'RSCRIPT'
options(repos = c(CRAN = "https://cloud.r-project.org"))
.libPaths(c("/home/srussell/Rlibs/4.5", .libPaths()))

required <- c("MOFA2", "UCell", "BayesPrism")
for (pkg in required) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat(sprintf("[%s] Installing %s\n", format(Sys.time()), pkg))
        if (pkg %in% c("MOFA2", "UCell")) {
            BiocManager::install(pkg, lib = "/home/srussell/Rlibs/4.5",
                                 update = FALSE, ask = FALSE)
        } else if (pkg == "BayesPrism") {
            # BayesPrism is on GitHub: Danko-Lab/BayesPrism
            if (!requireNamespace("remotes", quietly = TRUE)) {
                install.packages("remotes", lib = "/home/srussell/Rlibs/4.5")
            }
            remotes::install_github("Danko-Lab/BayesPrism/BayesPrism",
                                    lib = "/home/srussell/Rlibs/4.5",
                                    upgrade = "never")
        }
    } else {
        cat(sprintf("[%s] %s already present\n", format(Sys.time()), pkg))
    }
}

# Install Python mofapy2 backend for MOFA2 via R helper.
# This creates a per-user virtualenv under ~/.cache/R/basilisk or similar.
if (requireNamespace("MOFA2", quietly = TRUE)) {
    cat(sprintf("[%s] Setting up mofapy2 Python backend...\n", format(Sys.time())))
    tryCatch({
        MOFA2::install_mofapy2()
        cat("mofapy2 backend ready.\n")
    }, error = function(e) {
        cat("WARN: install_mofapy2() failed:", conditionMessage(e), "\n")
        cat("     Will fall back to manual virtualenv install in Stage G script.\n")
    })
}

# Final verification
need <- c("MOFA2", "UCell", "BayesPrism")
have <- need[need %in% rownames(installed.packages())]
miss <- setdiff(need, have)
cat("\n=== INSTALL VERIFICATION ===\n")
cat("have:", paste(have, collapse = " "), "\n")
cat("miss:", paste(miss, collapse = " "), "\n")
if (length(miss) > 0) quit(status = 1)
cat("All required packages installed.\n")
RSCRIPT

# Write manifest
MANIFEST="/project/def-srussell/nirosha_xeno/scorial/01_state/slurm_manifests/job_${SLURM_JOB_ID}.yaml"
cat > "$MANIFEST" <<EOF
---
job_id: ${SLURM_JOB_ID}
job_name: setup_install_libs
task_id: T-08-setup
stage: stage_G_setup
purpose: "Install MOFA2 + UCell + BayesPrism into ~/Rlibs/4.5; install mofapy2 Python backend via MOFA2::install_mofapy2()."
end_time: $(date -Iseconds)
node: $(hostname)
exit_code: 0
outputs:
  r_libs_path: /home/srussell/Rlibs/4.5
  packages_installed: [MOFA2, UCell, BayesPrism]
EOF
echo "Manifest written: $MANIFEST"
