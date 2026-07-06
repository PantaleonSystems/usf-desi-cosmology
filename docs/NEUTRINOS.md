# Tratamento de neutrinos massivos no pipeline USF

**Referência:** Tópico 2 do `pdfs/PLANO_CORRECAO_REFEREE1.md` (Comentário 1 do Referee 1).
**Decisão:** Σm_ν = 0.06 eV é incluído na contabilidade de densidades, com a
aproximação padrão de fluidos descrita abaixo, cujo erro foi **quantificado
contra o CAMB** (v1.6.6, m_ν = 0.06 eV, N_eff = 3.044) e é desprezível frente
às incertezas do DESI DR2.

## Convenções implementadas (`theory_module/background.py`)

1. **Baixo redshift (E(z), distâncias):** o parâmetro amostrado `Omega_m` é a
   fração de **matéria total** hoje (CDM + bárions + neutrinos
   não-relativísticos, ω_ν = Σm_ν/93.14 eV = 6.44×10⁻⁴). A radiação usa a
   densidade física fixa ω_r = ω_γ(1 + 0.2271·N_eff) com N_eff = 3.044
   (todas as espécies tratadas como relativísticas no termo (1+z)⁴).
2. **Alto redshift (integral de r_d):** para z > z_drag ≈ 1060 os neutrinos
   são ultrarrelativísticos (T_ν(z_drag) ≈ 0.18 eV ≫ m_ν/3 = 0.02 eV), então
   a integral do horizonte sonoro usa ω_cb = ω_m − ω_ν na matéria e o ω_r
   completo na radiação. A razão bárion/fóton usa **somente fótons**
   (ω_γ = 2.4728×10⁻⁵).

A dupla contagem residual (a espécie de 0.06 eV aparece como matéria em Ω_m e
como radiação em ω_r a baixo z) e a transição relativística→não-relativística
em z_nr ≈ 110, que a aproximação ignora, são as únicas fontes de erro.

## Erro medido contra o CAMB (parâmetros Planck 2018, m_ν = 0.06 eV)

| z | desvio em D_M(z) | desvio em H(z) |
|---|---|---|
| 0.10 | −0.0001% | +0.0002% |
| 0.51 | −0.0006% | +0.0014% |
| 1.00 | −0.0012% | +0.0027% |
| 1.49 | −0.0017% | +0.0040% |
| 2.33 | −0.0024% | +0.0059% |
| 5.00 | −0.0040% | +0.0112% |

Horizonte sonoro: r_d(pipeline) = 147.110 Mpc vs r_drag(CAMB) = 147.101 Mpc
→ **+0.006%**.

## Conclusão quantitativa

As medidas de BAO do DESI DR2 têm incertezas fracionárias de ~0.3–1% por
observável, e as distâncias do Pantheon+ ~1% por bin efetivo. O erro da
aproximação de neutrinos é ≤ 0.006% nos observáveis BAO (fator ≥ 50 abaixo
da estatística) e ≤ 0.004% nas distâncias de SNe. A aproximação é, portanto,
inteiramente adequada para esta análise; uma implementação com integrais
térmicas completas de Fermi–Dirac só se justificaria se Σm_ν fosse promovido
a parâmetro livre, o que está fora do escopo do artigo.

*Gerado na execução do Tópico 2 (2026-07-06); números reproduzíveis com
`tests/test_sound_horizon_vs_camb.py` e o snippet do log de execução em
`PROGRESSO_CODIGO.md`.*
