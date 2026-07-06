#!/usr/bin/env bash
set -euo pipefail

# USF pipeline runner — grade completa de runs do plano de correção
# (pdfs/PLANO_CORRECAO_REFEREE1.md, Tópicos 4/9/11).
#
# Ordem científica: os testes de validação SEMPRE rodam antes de amostrar
# (pré-condição do Tópico 9); as runs principais são as *_bao_sn (sem SH0ES);
# as variantes *_bao_sn_sh0es são rotuladas separadamente.

# Activate virtualenv if it exists
if [ -f ".venv/bin/activate" ]; then
    source .venv/bin/activate
fi

# ── Run grid (ordem: principais → variantes) ─────────────────────────────────
CONFIGS_MAIN=(usf_bao_sn lcdm_bao_sn)
CONFIGS_DATASETS=(usf_bao lcdm_bao usf_sn lcdm_sn)
CONFIGS_SH0ES=(usf_bao_sn_sh0es lcdm_bao_sn_sh0es)

# ── Parse flags ────────────────────────────────────────────────────────────────
CLEAN=false
SAMPLE=true
PLOT=true
TESTS=true
GRID="main"   # main | all | sh0es | datasets

usage() {
    echo "Usage: $0 [--clean] [--sample-only] [--plot-only] [--skip-tests] [--grid main|all|sh0es|datasets]"
    echo ""
    echo "  (no flags)       Run validation tests, main grid (BAO+SN), then plots"
    echo "  --grid all       Run the full grid (main + per-dataset + SH0ES variants)"
    echo "  --grid sh0es     Run only the *_bao_sn_sh0es variants"
    echo "  --grid datasets  Run only the per-dataset runs (BAO-only, SN-only)"
    echo "  --clean          Delete chains and results, then exit"
    echo "  --sample-only    Only run MCMC sampling (no plots)"
    echo "  --plot-only      Only generate plots (chains must already exist)"
    echo "  --skip-tests     Skip the pytest validation gate (NOT recommended)"
    exit 1
}

for arg in "$@"; do
    case $arg in
        --clean)       CLEAN=true; SAMPLE=false; PLOT=false; TESTS=false ;;
        --sample-only) PLOT=false ;;
        --plot-only)   SAMPLE=false; TESTS=false ;;
        --skip-tests)  TESTS=false ;;
        --grid)        ;; # value handled below
        main|all|sh0es|datasets) GRID=$arg ;;
        --help|-h)     usage ;;
        *) echo "Unknown flag: $arg"; usage ;;
    esac
done

case $GRID in
    main)     RUNS=("${CONFIGS_MAIN[@]}") ;;
    datasets) RUNS=("${CONFIGS_DATASETS[@]}") ;;
    sh0es)    RUNS=("${CONFIGS_SH0ES[@]}") ;;
    all)      RUNS=("${CONFIGS_MAIN[@]}" "${CONFIGS_DATASETS[@]}" "${CONFIGS_SH0ES[@]}") ;;
esac

# ── Step 0: Clean ──────────────────────────────────────────────────────────────
if [ "$CLEAN" = true ]; then
    echo "==> Cleaning chains and results..."
    for run in "${CONFIGS_MAIN[@]}" "${CONFIGS_DATASETS[@]}" "${CONFIGS_SH0ES[@]}"; do
        rm -f "chains/${run}."*
    done
    rm -f results/*.pdf results/*.png
    echo "Clean-only mode — exiting."
    exit 0
fi

# ── Step 1: Validation gate (Tópico 9: pré-condição para produção) ───────────
if [ "$TESTS" = true ]; then
    echo "==> Running validation suite (astropy/CAMB cross-checks)..."
    python -m pytest tests/ -q
    echo "    Validation passed."
fi

# ── Step 2: MCMC sampling ──────────────────────────────────────────────────────
if [ "$SAMPLE" = true ]; then
    for run in "${RUNS[@]}"; do
        echo ""
        echo "==> Running MCMC: cobaya_configs/${run}.yaml"
        cobaya-run "cobaya_configs/${run}.yaml" --resume
    done
    echo ""
    echo "    Sampling complete."
fi

# ── Step 3: Plots (a partir das runs principais por padrão) ───────────────────
if [ "$PLOT" = true ]; then
    USF_CHAIN="./chains/usf_bao_sn"
    LCDM_CHAIN="./chains/lcdm_bao_sn"
    if [ "$GRID" = "sh0es" ]; then
        USF_CHAIN="./chains/usf_bao_sn_sh0es"
        LCDM_CHAIN="./chains/lcdm_bao_sn_sh0es"
    fi

    echo ""
    echo "==> Generating figures from ${USF_CHAIN} / ${LCDM_CHAIN}..."
    python plot_scripts/plot_triangle.py    --usf-chain "$USF_CHAIN" --lcdm-chain "$LCDM_CHAIN"
    python plot_scripts/plot_desi_bao.py    --usf-chain "$USF_CHAIN" --lcdm-chain "$LCDM_CHAIN"
    python plot_scripts/plot_hubble_diag.py --usf-chain "$USF_CHAIN" --lcdm-chain "$LCDM_CHAIN"
    python plot_scripts/plot_omega_usf.py   --chain "$USF_CHAIN"

    echo ""
    echo "==> All figures saved to results/:"
    ls results/
fi

echo ""
echo "Pipeline finished."
