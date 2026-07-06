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

## Em Andamento

- [ ] Tópico 2 – Separar constantes de radiação em `usf_constants.py` (ω_γ vs ω_r)
- [ ] Tópico 2 – Usar ω_γ na razão R e ω_r no termo (1+z)⁴, nos dois módulos
- [ ] Tópico 2 – z_drag via Eisenstein & Hu (1998) em função de (ω_b, ω_m)
- [ ] Tópico 2 – Prior gaussiano de BBN em ω_b
- [ ] Tópico 2 – Tratamento de neutrinos massivos documentado

## Pendentes

- [ ] Tópico 3 – Auditar `sn.pantheonplus` (tratamento de M) e documentar
- [ ] Tópico 3 – Esquema canônico de calibração nos YAMLs (sem dupla contagem)
- [ ] Tópico 3 – Reescrever `plot_hubble_diag.py` (offset M marginalizado)
- [ ] Tópico 9 – Módulo único de background + refactor theory/likelihood/plots
- [ ] Tópico 9 – Validação cruzada vs CAMB no CI (tolerâncias documentadas)
- [ ] Tópico 9 – Testes de regressão de consistência interna
- [ ] Tópico 9 – Reprodutibilidade (pinagem, seeds, run_pipeline.sh)
- [ ] Tópico 9 – Renomear FEU → USF no código
- [ ] Tópico 5 – Módulo de CMB distance priors (R, ℓ_A, ω_b + covariância)
- [ ] Tópico 5 – Incluir CMB distance priors nos YAMLs de produção
- [ ] Tópico 6 – YAMLs de nested sampling (PolyChord)
- [ ] Tópico 6 – Script de geração de mocks e ajuste em lote
- [ ] Tópico 7 – Estender prior de z_trans nos YAMLs
- [ ] Tópico 11 – Atualizar README e `run_pipeline.sh`

## Log de execução

| Data | Tópico | Evento |
|---|---|---|
| 2026-07-06 | 1 | Suíte criada; 18 passed, 2 xfailed (bugs documentados). Desvio de r_d: **+4.401%** vs CAMB. χ² BAO no best-fit DESI: **227.5/13 pontos**. |
