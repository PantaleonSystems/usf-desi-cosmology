#!/usr/bin/env bash
set -euo pipefail

# Activate virtualenv if it exists
if [ -f ".venv/bin/activate" ]; then
    source .venv/bin/activate
fi

# ── Parse flags ────────────────────────────────────────────────────────────────
CLEAN=false
SAMPLE=true
PLOT=true

usage() {
    echo "Usage: $0 [--clean] [--sample-only] [--plot-only]"
    echo ""
    echo "  (no flags)     Clean chains+results, run MCMC, generate plots"
    echo "  --clean        Only delete chains and results, then exit"
    echo "  --sample-only  Only run MCMC sampling (no plots)"
    echo "  --plot-only    Only generate plots (chains must already exist)"
    exit 1
}

for arg in "$@"; do
    case $arg in
        --clean)       CLEAN=true; SAMPLE=false; PLOT=false ;;
        --sample-only) PLOT=false ;;
        --plot-only)   SAMPLE=false ;;
        --help|-h)     usage ;;
        *) echo "Unknown flag: $arg"; usage ;;
    esac
done

# ── Step 1: Clean ──────────────────────────────────────────────────────────────
echo "==> Cleaning previous chains and results..."
rm -f chains/feu_bao_sn_sh0es.* chains/lcdm_bao_sn_sh0es.*
rm -f results/*.pdf
echo "    Done."

if [ "$CLEAN" = true ]; then
    echo "Clean-only mode — exiting."
    exit 0
fi

# ── Step 2: MCMC sampling ──────────────────────────────────────────────────────
if [ "$SAMPLE" = true ]; then
    echo ""
    echo "==> Running FEU MCMC (cobaya_configs/feu_bao_sn_sh0es.yaml)..."
    cobaya-run cobaya_configs/feu_bao_sn_sh0es.yaml

    echo ""
    echo "==> Running ΛCDM MCMC (cobaya_configs/lcdm_bao_sn_sh0es.yaml)..."
    cobaya-run cobaya_configs/lcdm_bao_sn_sh0es.yaml

    echo ""
    echo "    Sampling complete."
fi

# ── Step 3: Plots ──────────────────────────────────────────────────────────────
if [ "$PLOT" = true ]; then
    echo ""
    echo "==> Generating triangle plot + H0 posterior..."
    python plot_scripts/plot_triangle.py

    echo ""
    echo "==> Generating BAO observables + comparison triangle..."
    python plot_scripts/plot_desi_bao.py

    echo ""
    echo "==> Generating Hubble diagram..."
    python plot_scripts/plot_hubble_diag.py

    echo ""
    echo "==> All figures saved to results/:"
    ls results/
fi

echo ""
echo "Pipeline finished."
