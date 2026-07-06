# Progresso da Execução das Tarefas [CÓDIGO]

**Plano:** `pdfs/PLANO_CORRECAO_REFEREE1.md`
**Última atualização:** 2026-07-06

## Tarefas Concluídas

### Tópico 1 — Auditoria e validação do pipeline

- [x] Tópico 1 – Criar suíte de testes `tests/` com pytest (E(z), D_C, D_A vs astropy `FlatLambdaCDM`, tolerância < 0.05%) → **18 testes PASSED** na grade (H0, Ω_m); limite ΛCDM (α_q=0) e normalização E(0)=1 verificados.
- [x] Tópico 1 – Testar cálculo de r_d contra CAMB (`tests/test_sound_horizon_vs_camb.py`) → **desvio documentado ANTES da correção: pipeline r_s = 153.575 Mpc vs CAMB rdrag = 147.101 Mpc (+4.401%)**. Teste marcado `xfail` até o Tópico 2.
- [x] Tópico 1 – Testar likelihood BAO isolada no best-fit publicado do DESI DR2 (Ω_m=0.2975, H0=68.5) → **χ² = 227.5 para 13 pontos** (esperado ~13). Teste marcado `xfail` até o Tópico 2.

**Diagnóstico da auditoria (explica o ΛCDM anômalo H0 = 64.2):**
a geometria ΛCDM do módulo está correta (< 0.05% vs astropy); o problema está
inteiramente no horizonte sonoro. O uso de `OMEGA_GAMMA_H2 = 4.15e-5` (radiação
**total**) como densidade de **fótons** na razão R = 3ω_b/(4ω_γ(1+z)) reduz R,
aumenta o sound speed e infla r_d em +4.4%. Pela degenerescência r_d·h das BAO,
H0 é puxado para baixo em ≈ 4.4%: 67.4 × (1 − 0.044) ≈ 64.4 ≈ **64.2** (o valor
anômalo da Tabela 2 do artigo). Este número vai para a carta-resposta.

### Tópico 2 — Setor de radiação e horizonte sonoro r_d

- [x] Tópico 2 – Separar constantes de radiação em `usf_constants.py`: `OMEGA_GAMMA_H2 = 2.4728e-5` (fótons) e `OMEGA_R_H2 = ω_γ(1+0.2271·3.044) ≈ 4.182e-5` (radiação total).
- [x] Tópico 2 – ω_γ usado na razão R do sound speed; ω_r no termo (1+z)⁴ — implementado uma única vez em `theory_module/background.py` (novo módulo compartilhado, antecipa o Tópico 9) e consumido por theory e likelihood.
- [x] Tópico 2 – z_drag(ω_b, ω_m) por ponto da cadeia. **Desvio do plano, documentado:** a fórmula EH98 dá z_d ≈ 1021 vs 1060 do CAMB (viés de +2.4% em r_d — trocaria um bug por outro). Implementada grade calibrada no CAMB (`theory_module/data/zdrag_camb_grid.npz`, gerada por `analysis_scripts/make_zdrag_grid.py`; também tabula z* para o Tópico 5) com interpolação; EH98 mantida como referência (`z_drag_eh98`) com teste documentando a discrepância. Cross-check independente com a fórmula de Aubourg et al. (2015) — concordância < 0.3%.
- [x] Tópico 2 – Prior gaussiano de BBN em ω_b (0.02218 ± 0.00055) nos dois YAMLs; ω_b agora é parâmetro amostrado (`omega_b`) da teoria e da likelihood BAO; r_d depende só de densidades físicas (teste `test_rd_independent_of_H0_and_usf`).
- [x] Tópico 2 – Neutrinos: Σm_ν = 0.06 eV em Ω_m a baixo z, relativísticos na integral de r_d (ω_cb = ω_m − ω_ν); erro quantificado contra CAMB em `docs/NEUTRINOS.md` (≤ 0.006% nos observáveis BAO).

**Resultado da correção:** r_d = 147.110 Mpc vs CAMB 147.101 (**+0.006%**, era +4.4%);
χ² BAO no best-fit DESI DR2 = **10.3 para 13 pontos** (era 227.5). Suíte: 27 passed.

### Tópico 3 — Calibração das supernovas (degenerescência M–H0)

- [x] Tópico 3 – Auditoria de `sn.pantheonplus` (Cobaya 3.6.2) documentada em `docs/SN_CALIBRATION.md`: usa `m_b_corr` (não calibrado), corte z > 0.01 (exclui calibradores Cefeidas) e **marginaliza M analiticamente** → a likelihood de SNe é insensível a H0; **não há dupla contagem** com o prior SH0ES no nível da likelihood.
- [x] Tópico 3 – Esquema canônico definido e documentado (`docs/SN_CALIBRATION.md`): runs principais `*_bao_sn` sem informação local de H0; variante rotulada `*_bao_sn_sh0es` com o prior gaussiano (nunca `pantheonplusshoes` + prior juntos; nunca `MU_SH0ES` em análises). Grade completa de YAMLs criada no Tópico 9.
- [x] Tópico 3 – `plot_hubble_diag.py` reescrito: `m_b_corr` + máscara da likelihood + offset M̂ ajustado por modelo (média ponderada, mimetiza a marginalização analítica); resíduos mostram apenas a *forma* de d_L(z). Smoke test: **com M marginalizado, ΛCDM e USF ajustam a forma das SNe igualmente bem (χ²/N diag 0.442 vs 0.446)** — confirma que o "drift" da Fig. 4 original era artefato da calibração SH0ES embutida em `MU_SH0ES`.

### Tópico 9 — Qualidade e unificação do código do pipeline

- [x] Tópico 9 – Módulo único `theory_module/background.py` (criado no Tópico 2) agora consumido por **todos**: theory (`usf_background.USF`), likelihood BAO, e os 4 plot scripts (nenhuma física local restante; `grep` confirma zero duplicação).
- [x] Tópico 9 – Likelihood BAO refatorada para o sistema **requirements/provider** do Cobaya (`get_requirements` → `provider.get_H/get_comoving_distance/get_r_s`); sem params próprios, sem física local.
- [x] Tópico 9 – Renomeado FEU → USF: `theory_module/usf_background.py` (classe `USF`), YAMLs, plots (`USF_H0_posterior.*`), labels. Zero referências residuais.
- [x] Tópico 9 – Grade completa de YAMLs (tarefa [CÓDIGO] do Tópico 4): `usf/lcdm × {bao_sn, bao, sn, bao_sn_sh0es}` com seeds fixos (1001–1008), SH0ES explícito (73.04 ± 1.04) apenas nas variantes `_sh0es`.
- [x] Tópico 9 – Testes de regressão (`test_regression_consistency.py`): theory ≡ background bit a bit; limites Ω_USF; **teste de integração Cobaya** (`test_cobaya_integration.py`): modelo completo montado do YAML, posterior finito, χ²_BAO via provider < 40, limite aninhado (α_q=0 independe de z_trans no modelo completo).
- [x] Tópico 9 – Validação cruzada externa (`test_camb_crossval.py`): D_M, D_A, H(z) < 0.1% e r_d < 0.2% vs CAMB em 4 cosmologias cobrindo o volume de prior — gate obrigatório no `run_pipeline.sh` antes de amostrar.
- [x] Tópico 9 – Reprodutibilidade: `requirements.txt` com astropy/camb/pytest pinados; seeds nos YAMLs; `run_pipeline.sh` reescrito (gate de testes + grade `--grid main|all|sh0es|datasets` + plots parametrizados).

### Tópico 5 — Consistência com o CMB (distance priors)

- [x] Tópico 5 – Módulo `theory_module/cmb_distance_priors.py`: vetor (R, ℓ_A, ω_b) = (1.7502±0.0046, 301.471±0.090, 0.02236±0.00015) com matriz de correlação completa de **Chen, Huang & Wang 2019 (arXiv:1808.05724, Tabela I, verificada no PDF)**. Mapeamento teoria→dado reproduz **exatamente** a receita de compressão deles (Eqs. 1–6: z* das cadeias Planck/CAMB via nossa grade; E(z) com Ω_r = Ω_m/(1+z_eq); r_s com 31500ω_b(T/2.7)⁻⁴) — validado no mean do Planck 2018: predição (1.7503, 301.429) vs dado (1.7502, 301.471), **χ² = 0.29**. Termo USF entra na linha de visada (teste garante sensibilidade a α_q — exatamente o efeito que o revisor exigiu constranger).
- [x] Tópico 5 – Bônus de precisão: `bg.comoving_distance_cmb` com transição de neutrino massivo resolvida (fitting de Komatsu+) → D_M(z*) a **−0.0007%** do CAMB (era −0.027% na aproximação de fluido); grade de z* já tabulada junto com z_drag; provider `get_rs_zstar`/`get_z_star`/`get_DM` expostos na teoria.
- [x] Tópico 5 – 4 YAMLs de produção `+cmb` criados (`usf/lcdm × {bao_sn, bao_sn_sh0es} _cmb`, seeds 1009–1012) e grupo `--grid cmb` no `run_pipeline.sh` (incluído em `--grid all`). Suíte: 44 passed.
### Tópico 6 — Estatística de comparação de modelos

- [x] Tópico 6 – 4 YAMLs de nested sampling (`*_polychord.yaml`: usf/lcdm × {bao_sn, bao_sn_cmb}, nlive=500, precision 0.001, seeds 2001–2004). Requer `cobaya-install polychord` antes das runs [ANÁLISE].
- [x] Tópico 6 – `analysis_scripts/run_mocks.py`: calibração frequentista do Δχ² sob a nula ΛCDM — mocks de BAO (cov DESI) + SNe (cov STAT+SYS completa, M marginalizado por projeção idêntica à do Cobaya), ajuste de ΛCDM e USF **perfilando z_trans em grade geomspace(0.5, 20)** (nunca fixo — look-elsewhere), prior BBN incluso, flag `--with-sh0es` (sorteia o mock do prior local). Smoke test 2 mocks: Δχ² = {1.06, 3.76} sob a nula — mecânica validada. Produção: `--n-mocks 300+` ([ANÁLISE]).

### Tópico 7 — Prior de z_trans

- [ ] Tópico 7 – Estender prior de z_trans nos YAMLs

### Tópico 11 — Empacotamento

- [ ] Tópico 11 – Atualizar README

## Log de execução

| Data | Tópico | Evento |
|---|---|---|
| 2026-07-06 | 1 | Suíte criada; 18 passed, 2 xfailed (bugs documentados). Desvio de r_d: **+4.401%** vs CAMB. χ² BAO no best-fit DESI: **227.5/13 pontos**. |
| 2026-07-06 | 2 | Setor de radiação corrigido; z_drag calibrado em grade CAMB; prior BBN em ω_b; r_d **+0.006%** vs CAMB; χ² BAO **10.3/13**. 27 passed. |
| 2026-07-06 | 3 | Auditoria sn.pantheonplus (M marginalizado analiticamente; sem dupla contagem); plot_hubble_diag.py refeito com M̂ por modelo — drift da Fig. 4 confirmado como artefato. |
| 2026-07-06 | 9 | Unificação completa: provider Cobaya, FEU→USF, grade de 8 YAMLs com seeds, cross-val CAMB como gate, 39 passed. |
| 2026-07-06 | 5 | CMB distance priors (Chen+19 verificado no paper) implementados com receita de compressão exata; χ²=0.29 no mean do Planck; 4 YAMLs +cmb. 44 passed. |
| 2026-07-06 | 6 | PolyChord YAMLs (4) + run_mocks.py (perfil em z_trans, cov SN completa, M projetado); smoke 2 mocks OK. |
