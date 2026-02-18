#!/bin/bash
# SLURM submission wrapper for merge-multisample-vcf.
# Auto-detects BIH HPC, Charité HPC, or local execution.
#
# Usage:
#   sbatch scripts/run_snakemake.sh                        # default config
#   sbatch scripts/run_snakemake.sh config/my_config.yaml  # custom config
#   bash   scripts/run_snakemake.sh --profile profiles/local --dry-run
#
# Extra arguments are forwarded verbatim to snakemake (e.g. --dry-run, --notemp).
#SBATCH --job-name=merge_multisample_vcf
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=2-00:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --output=slurm_logs/%x-%j.log

set -euo pipefail

# ---------------------------------------------------------------------------
# Cluster auto-detection
# ---------------------------------------------------------------------------
detect_cluster() {
    local fqdn
    fqdn=$(hostname -f 2>/dev/null || hostname)
    if [[ -d "/etc/xdg/snakemake/cubi-v1" ]] || [[ "$fqdn" =~ cubi|bihealth ]]; then
        echo "bih"
    elif [[ -f "/etc/profile.d/conda.sh" ]] || [[ "$fqdn" =~ charite|\.sc- ]]; then
        echo "charite"
    else
        echo "local"
    fi
}

CLUSTER=$(detect_cluster)
CONFIG_FILE="${1:-config/config.yaml}"
shift 1 2>/dev/null || true

# ---------------------------------------------------------------------------
# Conda activation (Charité requires explicit sourcing)
# ---------------------------------------------------------------------------
if [[ "$CLUSTER" == "charite" ]] && [[ -f /etc/profile.d/conda.sh ]]; then
    source /etc/profile.d/conda.sh
fi
conda activate snakemake 2>/dev/null || true

# ---------------------------------------------------------------------------
# TMPDIR — use cluster scratch where available
# ---------------------------------------------------------------------------
if [[ "$CLUSTER" == "bih" ]]; then
    BASE_TMPDIR="${HOME}/scratch/tmp"
else
    BASE_TMPDIR="${TMPDIR:-/tmp}/snakemake"
fi
mkdir -p "${BASE_TMPDIR}"
TMPDIR=$(mktemp -d "${BASE_TMPDIR}/sm.XXXXXX")
export TMPDIR
trap 'rm -rf "${TMPDIR}"' EXIT

mkdir -p slurm_logs

# ---------------------------------------------------------------------------
# Launch
# ---------------------------------------------------------------------------
echo "=== merge-multisample-vcf ==="
echo "  Cluster:  ${CLUSTER}"
echo "  Config:   ${CONFIG_FILE}"
echo "  Profile:  profiles/${CLUSTER}"
echo "  TMPDIR:   ${TMPDIR}"
echo "  Extra:    $*"
echo "  Start:    $(date)"
echo "============================="

snakemake \
    -s workflow/Snakefile \
    --configfile "${CONFIG_FILE}" \
    --workflow-profile profiles/default \
    --profile "profiles/${CLUSTER}" \
    "$@"

echo "=== Finished: $(date) ==="
