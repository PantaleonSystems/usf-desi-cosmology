#!/usr/bin/env bash
set -euo pipefail

# USF pipeline runner — grade completa de runs do plano de correção
# (pdfs/PLANO_CORRECAO_REFEREE1.md, Tópicos 4/9/11).
#
# Ordem científica: os testes de validação SEMPRE rodam antes de amostrar
# (pré-condição do Tópico 9); as runs principais são as *_bao_sn (sem SH0ES);
# as variantes *_bao_sn_sh0es são rotuladas separadamente.
#
# Layout de saídas:
#   chains/<run>.*              amostras + metadados do Cobaya
#   logs/<run>.log              log de cada run (um arquivo por run)
#   results/figures/<tag>/      figuras, um subdiretório por par de cadeias
#   results/data/               produtos numéricos (mocks, tabelas) — NÃO
#                               apagados por --clean, só por --clean-all

# Activate virtualenv if it exists
if [ -f ".venv/bin/activate" ]; then
    source .venv/bin/activate
fi

LOG_DIR="logs"
FIG_DIR="results/figures"
DATA_DIR="results/data"

# ── Run grid (ordem: principais → variantes) ─────────────────────────────────
CONFIGS_MAIN=(usf_bao_sn lcdm_bao_sn)
CONFIGS_DATASETS=(usf_bao lcdm_bao usf_sn lcdm_sn)
CONFIGS_SH0ES=(usf_bao_sn_sh0es lcdm_bao_sn_sh0es)
CONFIGS_CMB=(usf_bao_sn_cmb lcdm_bao_sn_cmb usf_bao_sn_sh0es_cmb lcdm_bao_sn_sh0es_cmb)
CONFIGS_POLYCHORD=(usf_bao_sn_polychord lcdm_bao_sn_polychord
                   usf_bao_sn_cmb_polychord lcdm_bao_sn_cmb_polychord)

# Todas as runs conhecidas — usado por --clean e pela checagem de nomes.
ALL_CONFIGS=("${CONFIGS_MAIN[@]}" "${CONFIGS_DATASETS[@]}" "${CONFIGS_SH0ES[@]}"
             "${CONFIGS_CMB[@]}" "${CONFIGS_POLYCHORD[@]}")

# ── Parse flags ────────────────────────────────────────────────────────────────
CLEAN=false
CLEAN_ALL=false
SAMPLE=true
PLOT=true
TESTS=true
FRESH=false
GRID="main"        # main | all | sh0es | datasets | cmb | polychord
CUSTOM_RUNS=""     # --runs a,b,c  (sobrepõe --grid; aceita YAMLs novos)
FIG_TAG=""         # --tag nome    (subdiretório em results/figures/)

usage() {
    cat <<'USAGE'
Usage: ./run_pipeline.sh [opções]

  (sem flags)      Testes de validação → grade principal (BAO+SN) → figuras

Seleção de runs:
  --grid main|all|sh0es|datasets|cmb|polychord
                   Grade pré-definida (default: main)
  --runs a,b,c     Lista explícita de runs (nomes de cobaya_configs/<nome>.yaml).
                   Use isto para testar parametrizações novas sem tocar nas grades.
  --tag NOME       Subdiretório das figuras: results/figures/NOME
                   (default: o nome da grade, ou 'custom' com --runs)

Execução:
  --fresh          Apaga chains/<run>.* antes de amostrar, em vez de --resume.
                   Necessário sempre que o YAML mudar de teoria/likelihood —
                   é o que evita o erro "Old and new run information not
                   compatible! Resuming not possible!".
  --sample-only    Só amostra (sem figuras)
  --plot-only      Só gera figuras (cadeias já devem existir)
  --skip-tests     Pula o gate de validação pytest (NÃO recomendado)

Limpeza:
  --clean          Apaga chains/, logs/ e results/figures/ e sai.
                   Preserva results/data/ (mocks são caros de refazer).
  --clean-all      Idem, incluindo results/data/.
USAGE
    exit 1
}

while [ $# -gt 0 ]; do
    case $1 in
        --clean)       CLEAN=true; SAMPLE=false; PLOT=false; TESTS=false ;;
        --clean-all)   CLEAN=true; CLEAN_ALL=true; SAMPLE=false; PLOT=false; TESTS=false ;;
        --sample-only) PLOT=false ;;
        --plot-only)   SAMPLE=false; TESTS=false ;;
        --skip-tests)  TESTS=false ;;
        --fresh)       FRESH=true ;;
        --grid)        shift; GRID="${1:-}" ;;
        --runs)        shift; CUSTOM_RUNS="${1:-}" ;;
        --tag)         shift; FIG_TAG="${1:-}" ;;
        # forma legada: `--grid all` já consumida acima; aceita `all` solto
        main|all|sh0es|datasets|cmb|polychord) GRID=$1 ;;
        --help|-h)     usage ;;
        *) echo "Unknown flag: $1"; usage ;;
    esac
    shift
done

if [ -n "$CUSTOM_RUNS" ]; then
    IFS=',' read -r -a RUNS <<< "$CUSTOM_RUNS"
    [ -z "$FIG_TAG" ] && FIG_TAG="custom"
else
    case $GRID in
        main)      RUNS=("${CONFIGS_MAIN[@]}") ;;
        datasets)  RUNS=("${CONFIGS_DATASETS[@]}") ;;
        sh0es)     RUNS=("${CONFIGS_SH0ES[@]}") ;;
        cmb)       RUNS=("${CONFIGS_CMB[@]}") ;;
        polychord) RUNS=("${CONFIGS_POLYCHORD[@]}") ;;
        all)       RUNS=("${CONFIGS_MAIN[@]}" "${CONFIGS_DATASETS[@]}"
                         "${CONFIGS_SH0ES[@]}" "${CONFIGS_CMB[@]}") ;;
        *) echo "Grade desconhecida: $GRID"; usage ;;
    esac
    [ -z "$FIG_TAG" ] && FIG_TAG="$GRID"
fi

OUTDIR="${FIG_DIR}/${FIG_TAG}"

# ── Step 0: Clean ──────────────────────────────────────────────────────────────
if [ "$CLEAN" = true ]; then
    echo "==> Limpando cadeias, logs e figuras..."
    for run in "${ALL_CONFIGS[@]}"; do
        rm -f "chains/${run}."*
        rm -rf "chains/${run}"          # PolyChord grava um diretório
    done
    rm -rf "$LOG_DIR" "$FIG_DIR"
    if [ "$CLEAN_ALL" = true ]; then
        echo "==> Limpando também ${DATA_DIR}/ (mocks, tabelas)..."
        rm -rf "$DATA_DIR"
    else
        echo "    ${DATA_DIR}/ preservado (use --clean-all para apagar também)."
    fi
    echo "Clean-only mode — exiting."
    exit 0
fi

# ── Checagem: todo run pedido precisa de um YAML ──────────────────────────────
for run in "${RUNS[@]}"; do
    if [ ! -f "cobaya_configs/${run}.yaml" ]; then
        echo "ERRO: cobaya_configs/${run}.yaml não existe."
        exit 1
    fi
done

mkdir -p "$LOG_DIR" "$DATA_DIR"

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
        echo "==> Running MCMC: cobaya_configs/${run}.yaml  (log: ${LOG_DIR}/${run}.log)"
        if [ "$FRESH" = true ]; then
            rm -f "chains/${run}."*
            rm -rf "chains/${run}"
            cobaya-run "cobaya_configs/${run}.yaml" 2>&1 | tee "${LOG_DIR}/${run}.log"
        else
            cobaya-run "cobaya_configs/${run}.yaml" --resume 2>&1 | tee "${LOG_DIR}/${run}.log"
        fi
    done
    echo ""
    echo "    Sampling complete."
fi

# ── Step 3: Plots (a partir das runs principais por padrão) ───────────────────
if [ "$PLOT" = true ]; then
    USF_CHAIN="./chains/usf_bao_sn"
    LCDM_CHAIN="./chains/lcdm_bao_sn"
    case $GRID in
        sh0es) USF_CHAIN="./chains/usf_bao_sn_sh0es"; LCDM_CHAIN="./chains/lcdm_bao_sn_sh0es" ;;
        cmb)   USF_CHAIN="./chains/usf_bao_sn_cmb";   LCDM_CHAIN="./chains/lcdm_bao_sn_cmb" ;;
    esac
    # Com --runs, usa o primeiro par usf/lcdm da lista.
    if [ -n "$CUSTOM_RUNS" ]; then
        for run in "${RUNS[@]}"; do
            case $run in
                usf_*)  [ "$USF_CHAIN" = "./chains/usf_bao_sn" ] && USF_CHAIN="./chains/${run}" ;;
                lcdm_*) [ "$LCDM_CHAIN" = "./chains/lcdm_bao_sn" ] && LCDM_CHAIN="./chains/${run}" ;;
            esac
        done
    fi

    mkdir -p "$OUTDIR"
    echo ""
    echo "==> Generating figures from ${USF_CHAIN} / ${LCDM_CHAIN} into ${OUTDIR}/..."
    python plot_scripts/plot_triangle.py    --usf-chain "$USF_CHAIN" --lcdm-chain "$LCDM_CHAIN" --outdir "$OUTDIR"
    python plot_scripts/plot_desi_bao.py    --usf-chain "$USF_CHAIN" --lcdm-chain "$LCDM_CHAIN" --outdir "$OUTDIR"
    python plot_scripts/plot_hubble_diag.py --usf-chain "$USF_CHAIN" --lcdm-chain "$LCDM_CHAIN" --outdir "$OUTDIR"
    python plot_scripts/plot_omega_usf.py   --chain "$USF_CHAIN" --outdir "$OUTDIR"

    echo ""
    echo "==> All figures saved to ${OUTDIR}/:"
    ls "$OUTDIR"
fi

echo ""
echo "Pipeline finished."
