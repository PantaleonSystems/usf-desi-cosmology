# Calibração das supernovas: auditoria e esquema canônico

**Referência:** Tópico 3 do `pdfs/PLANO_CORRECAO_REFEREE1.md` (Comentário 3 do Referee 1).

## O que a likelihood `sn.pantheonplus` do Cobaya faz com M (auditoria, Cobaya 3.6.2)

Fonte auditada: `cobaya/likelihoods/sn/pantheonplus.py`.

1. Lê **`m_b_corr`** (magnitudes aparentes corrigidas, **não calibradas** pela
   SH0ES), `zHD` e `zHEL` do arquivo `Pantheon+SH0ES.dat`, com a covariância
   completa STAT+SYS.
2. Aplica o corte **z_HD > 0.01**, que remove as SNe hospedadas em galáxias
   com calibração de Cefeidas (os "calibrators" da SH0ES ficam de fora).
3. Com `use_abs_mag: False` (default, e o que usamos), a magnitude absoluta M
   é **marginalizada analiticamente**: o método `_marginalize_abs_mag()`
   projeta o modo constante para fora da matriz de covariância inversa
   (equivalente à marginalização gaussiana exata sobre M com prior plano).

**Consequência direta (degenerescência M–H0):** com M marginalizado, a
likelihood de SNe é **insensível a H0** — ela restringe apenas a *forma* de
d_L(z), i.e. Ω_m, α_q, z_trans. Portanto:

- ✅ **Não há dupla contagem** entre `sn.pantheonplus` e o prior SH0ES em H0
  no nível da likelihood: toda a informação local de H0 entra uma única vez,
  pelo prior gaussiano.
- ❌ A narrativa da Sec. 4.4 do artigo (drift dos resíduos ΛCDM "causado pelo
  H0 baixo") é insustentável: qualquer offset global entre modelo e dados é
  absorvido por M. O drift da Fig. 4 original era artefato de comparar
  `MU_SH0ES` (μ **calibrado pela SH0ES**, coluna 11 do catálogo) com modelos
  sem nenhum ajuste de offset — exatamente a circularidade apontada pelo
  revisor. Corrigido em `plot_scripts/plot_hubble_diag.py`.

Para referência: `sn.pantheonplusshoes` (não usada por nós) inclui os
calibradores e usa as distâncias de Cefeidas como teoria para eles; nesse
caso M (`Mb`) precisa ser amostrado (`use_abs_mag: True`) e um prior SH0ES em
H0 **não pode** ser combinado (dupla contagem).

## Esquema canônico adotado (YAMLs em `cobaya_configs/`)

| Rótulo | Likelihoods | Papel de M | Informação local de H0 |
|---|---|---|---|
| `*_bao_sn` | DESI BAO + `sn.pantheonplus` (+ prior BBN em ω_b) | marginalizado analiticamente | **nenhuma** (runs principais) |
| `*_bao_sn_sh0es` | os anteriores + prior gaussiano H0 = 73.04 ± 1.04 | marginalizado analiticamente | via prior em H0, uma única vez |

Regras invioláveis:

1. Nunca combinar `sn.pantheonplusshoes` com o prior SH0ES em H0.
2. Nunca usar a coluna `MU_SH0ES` do catálogo em análises ou figuras de
   comparação de modelos — usar `m_b_corr` com offset M ajustado/
   marginalizado por modelo.
3. O prior SH0ES entra apenas em runs explicitamente rotulados `_sh0es`,
   apresentados como variante, nunca como resultado principal.

## Tratamento de M nas figuras

`plot_scripts/plot_hubble_diag.py` reescrito: usa `m_b_corr` com o corte
z > 0.01 da likelihood e, para cada modelo, ajusta o offset global
M̂ = média ponderada de (m_b_corr − μ_model) (pesos 1/σ², erros diagonais —
suficiente para visualização; a inferência usa a covariância completa via
Cobaya). Os resíduos mostrados são Δ = m_b_corr − μ_model − M̂, que refletem
apenas a *forma* de d_L(z) — comparável entre modelos e imune à
degenerescência M–H0.
