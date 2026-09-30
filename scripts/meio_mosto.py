"""
meio_mosto.py -- camada de MEIO: mosto sintetico MS300 -> bounds de captacao.

Escrito de raiz a 2026-09-25. Nao herda de `mosto_medium.py` nem de
`ms300_medium_fixed.py`; a receita e a dissociacao vem de la porque sao DADOS
auditados, e estao creditadas onde aparecem.

A equacao
---------
Relacao de Sauer para cultura em batch. Em batch, dS/dt = -qS.X e dX/dt = mu.X,
logo dS/dX = -qS/mu e ΔS = (qS/mu).ΔX. Com X(t1) << X(t2):

    qS = mu . ΔS / X(t2)                      [mmol/gDW/h]

    ΔS     concentracao disponivel no meio     [mmol/L]
    mu     taxa especifica de crescimento      [1/h]
    X(t2)  biomassa no fim do crescimento      [gDW/L]

O tempo nao entra. Quem ocupa o lugar de `t` sao mu e a biomassa final, que sao
grandezas medidas. A versao anterior usava qS = S0/(X.t) com X = 0.5 gDW/L e
t = 24 h escolhidos sem fonte, o que da um denominador de 12 gDW.h/L contra os
29 que uma fermentacao vinica real tem.

Verificacao que se deve fazer a qualquer bound deste modulo: multiplicar de
volta. qS . X(t2) / mu tem de dar a concentracao da receita. E o que
`verificar_bound()` faz.

Os dois parametros
------------------
mu = 0.2 h-1 e X(t2) = 5.8 gDW/L, ambos de Varela, Pizarro & Agosin (2004),
Appl Environ Microbiol 70:3392-3400: mosto sintetico, 240 g/L de acucar,
300 mg N/L de azoto assimilavel, que e o azoto do MS300.

Um so par para os quatro modelos. Em qS = mu.ΔS/X(t2) o par descreve o REGIME
da fermentacao, nao o organismo: o bound significa "o que esta fermentacao
disponibiliza por grama de celulas por hora". Se cada modelo usasse o seu
proprio mu, o bound dependia da previsao que ele serve para constranger.

Sobre as fracoes de eficiencia
------------------------------
O modulo `legume_medium1_v2.py` do projecto anterior tinha uma fraccao por
composto (acucares 1.0, proteina 0.20, lipidos 0.12) porque numa matriz solida
a proteina esta ligada e so parte fica acessivel. No mosto sintetico esta tudo
dissolvido e a fraccao e 1.0 por construcao. O parametro fica no sitio, com o
valor 1.0 e esta nota, para a estrutura ser a mesma e para quem quiser testar
uma hipotese de acessibilidade ter onde a por.
"""

from __future__ import annotations

import pandas as pd

# ═══════════════════════════════════════════════════════════ 1 · parametros ══
MU_REF = 0.2        # h-1     Varela, Pizarro & Agosin 2004, AEM 70:3392
X_REF = 5.8         # gDW/L   idem, biomassa final a 300 mg N/L

# Cenarios: variam o par de parametros, nao a composicao. O mosto sintetico tem
# composicao fixa por definicao; o que tem incerteza e o regime de fermentacao.
#   minimo  -- fermentacao lenta, 50 mg N/L, X = 1.5 gDW/L (Varela et al. 2004)
#   central -- fermentacao normal, 300 mg N/L (Varela et al. 2004)
#   maximo  -- MS300 a 324 mg N/L, X proximo de 8 g/L (Brou et al., OENO One)
CENARIOS = {
    "minimo":  {"mu": 0.05, "x_final": 1.5},
    "central": {"mu": 0.20, "x_final": 5.8},
    "maximo":  {"mu": 0.25, "x_final": 8.0},
}

EFICIENCIA_OMISSAO = 1.0   # meio liquido: tudo dissolvido (ver docstring)

# Coeficiente de `ash 1 g` na BIOMASS dos tres modelos CarveFungi.
ASH_POR_GDW = 0.048


# ═════════════════════════════════════════════════════════════ 2 · receita ══
# Evers et al. 2023, Foods 12(5):972, seccao 2.2.
# MW e RECIPE_MG_L transcritos de `ms300_medium_fixed.py`, que os auditou.
MW = {
    "glucose": 180.16, "fructose": 180.16, "DL_malate": 134.09, "citrate": 192.12,
    "KH2PO4": 136.09, "K2SO4": 174.26, "MgSO4_7H2O": 246.47, "CaCl2_2H2O": 147.01,
    "NaCl": 58.44, "NH4Cl": 53.49,
    "MnSO4_H2O": 169.02, "ZnSO4_7H2O": 287.56, "CuSO4_5H2O": 249.69, "KI": 166.00,
    "CoCl2_6H2O": 237.93, "H3BO3": 61.83, "NH4_6Mo7O24": 1163.9,
    "ergosterol": 396.65, "oleate": 282.46, "tween80": 1310.0,
    "myo_inositol": 180.16, "pantothenate_hemiCa": 476.53, "nicotinate": 123.11,
    "pyridoxine": 169.18, "thiamine_HCl": 337.27, "biotin": 244.31,
    "tyrosine": 181.19, "tryptophan": 204.23, "isoleucine": 131.17,
    "aspartate": 133.10, "glutamate": 147.13, "arginine": 174.20,
    "leucine": 131.17, "threonine": 119.12, "glycine": 75.07,
    "glutamine": 146.14, "alanine": 89.09, "valine": 117.15,
    "methionine": 149.21, "phenylalanine": 165.19, "serine": 105.09,
    "histidine": 155.15, "lysine": 146.19, "cysteine": 121.16,
}

RECEITA_MG_L = {
    "glucose": 100000.0, "fructose": 100000.0, "DL_malate": 6000.0, "citrate": 6000.0,
    "KH2PO4": 750.0, "K2SO4": 500.0, "MgSO4_7H2O": 250.0, "CaCl2_2H2O": 155.0,
    "NaCl": 200.0, "NH4Cl": 460.0,
    "MnSO4_H2O": 4.0, "ZnSO4_7H2O": 4.0, "CuSO4_5H2O": 1.0, "KI": 1.0,
    "CoCl2_6H2O": 0.4, "H3BO3": 1.0, "NH4_6Mo7O24": 1.0,
    "ergosterol": 3.0, "oleate": 1.0, "tween80": 100.0,
    "myo_inositol": 20.0, "pantothenate_hemiCa": 1.5, "nicotinate": 2.0,
    "pyridoxine": 0.25, "thiamine_HCl": 0.25, "biotin": 0.003,
    "tyrosine": 18.0, "tryptophan": 179.0, "isoleucine": 33.0,
    "aspartate": 45.0, "glutamate": 120.0, "arginine": 374.0,
    "leucine": 48.0, "threonine": 76.0, "glycine": 18.0,
    "glutamine": 505.0, "alanine": 145.0, "valine": 45.0,
    "methionine": 31.0, "phenylalanine": 38.0, "serine": 79.0,
    "histidine": 33.0, "lysine": 17.0, "cysteine": 13.0,
    # Prolina ausente de proposito: nao consta da seccao 2.2 do Evers, e em
    # anaerobiose nao e assimilavel (prolina oxidase mitocondrial, depende de O2).
}

# Sais e Tween dissociam-se antes de virar bound. O Tween 80 e monooleato de
# polioxietileno(20)-sorbitano: 1 mol de oleato por mol nominal. Sem isto o
# oleato do MS300 fica 22x subestimado e a corrida para por falta de acido
# gordo insaturado, nao por falta de azoto.
DISSOCIACAO = {
    "KH2PO4": {"potassium": 1, "phosphate": 1},
    "K2SO4": {"potassium": 2, "sulfate": 1},
    "MgSO4_7H2O": {"magnesium": 1, "sulfate": 1},
    "CaCl2_2H2O": {"calcium": 1, "chloride": 2},
    "NaCl": {"sodium": 1, "chloride": 1},
    "NH4Cl": {"ammonium": 1, "chloride": 1},
    "MnSO4_H2O": {"manganese": 1, "sulfate": 1},
    "ZnSO4_7H2O": {"zinc": 1, "sulfate": 1},
    "CuSO4_5H2O": {"copper": 1, "sulfate": 1},
    "KI": {"potassium": 1, "iodide": 1},
    "CoCl2_6H2O": {"cobalt": 1, "chloride": 2},
    "H3BO3": {"borate": 1},
    "NH4_6Mo7O24": {"ammonium": 6, "molybdate": 7},
    "pantothenate_hemiCa": {"pantothenate": 2, "calcium": 1},
    "thiamine_HCl": {"thiamine": 1, "chloride": 1},
    "tween80": {"oleate": 1},
    "DL_malate": {"malate": 1},
}

# Atomos de azoto por molecula. A convencao enologica do YAN conta 1 N por
# aminoacido (alfa-amino) e da os 301 mg N/L que baptizam o MS300. A biomassa
# consome ATOMOS, e a arginina tem 4.
N_POR_MOLECULA = {
    "ammonium": 1, "arginine": 4, "histidine": 3, "tryptophan": 2,
    "glutamine": 2, "asparagine": 2, "lysine": 2,
    "alanine": 1, "aspartate": 1, "cysteine": 1, "glutamate": 1, "glycine": 1,
    "isoleucine": 1, "leucine": 1, "methionine": 1, "phenylalanine": 1,
    "proline": 1, "serine": 1, "threonine": 1, "tyrosine": 1, "valine": 1,
}
AMINOACIDOS = tuple(k for k in N_POR_MOLECULA if k != "ammonium")


def composicao_mmol_L(receita=None):
    """Receita em mg/L -> concentracoes em mmol/L, com os sais dissociados."""
    receita = RECEITA_MG_L if receita is None else receita
    out: dict[str, float] = {}
    for composto, mg_L in receita.items():
        mmol = mg_L / MW[composto]
        for especie, n in DISSOCIACAO.get(composto, {composto: 1}).items():
            out[especie] = out.get(especie, 0.0) + n * mmol
    return out


def yan(comp=None):
    """YAN nas duas convencoes, em mg N/L. A enologica e a que da o '300'."""
    comp = composicao_mmol_L() if comp is None else comp
    alfa = sum(v for k, v in comp.items() if k in N_POR_MOLECULA) * 14.007
    atomos = sum(v * N_POR_MOLECULA[k] for k, v in comp.items()
                 if k in N_POR_MOLECULA) * 14.007
    return {"YAN_alfa_amino_mgN_L": round(alfa, 1),
            "N_em_atomos_mgN_L": round(atomos, 1)}


# ═══════════════════════════════════════════════════════ 3 · concentracao → bound ══
def calcular_bound(mmol_L, mu=MU_REF, x_final=X_REF, eficiencia=EFICIENCIA_OMISSAO):
    """qS = mu . ΔS / X(t2), devolvido negativo (captacao), em mmol/gDW/h."""
    if mmol_L <= 0 or eficiencia <= 0:
        return 0.0
    return -(mu * mmol_L * eficiencia) / x_final


def verificar_bound(bound, mu=MU_REF, x_final=X_REF):
    """Multiplica de volta: |qS| . X(t2) / mu tem de dar a concentracao do meio.

    E a verificacao que responde a 'estes bounds nao sao pequenos de mais?'.
    O bound nao e pequeno nem grande: e a taxa que esvazia o meio exactamente
    quando a biomassa chega a X(t2)."""
    return abs(bound) * x_final / mu


def bounds_do_meio(comp=None, mu=MU_REF, x_final=X_REF, eficiencia=None):
    """{especie: bound negativo} para toda a composicao."""
    comp = composicao_mmol_L() if comp is None else comp
    eficiencia = eficiencia or {}
    return {k: calcular_bound(v, mu, x_final,
                              eficiencia.get(k, EFICIENCIA_OMISSAO))
            for k, v in comp.items()}


# ═════════════════════════════════════════════════════ 4 · mapa de trocas ══
# Sinonimos aceites para o NOME do metabolito da troca. As duas familias de
# ficheiros nomeiam de maneira diferente: o yeast-GEM escreve "L-glycine" e
# "L-arginine", o CarveFungi escreve "glycine zwitterion" e "L-argininium".
# E esta tabela, e nao um mapa por modelo, o artefacto a auditar.
SINONIMOS = {
    "glucose":       ["D-glucose", "alpha-D-glucose"],
    "fructose":      ["D-fructose"],
    "malate":        ["(S)-malate", "L-malate", "malate"],
    "citrate":       ["citrate", "citrate(3-)"],
    "potassium":     ["potassium", "K(+)"],
    "phosphate":     ["phosphate", "hydrogenphosphate"],
    "sulfate":       ["sulphate", "sulfate"],
    "magnesium":     ["Mg(2+)", "magnesium"],
    "calcium":       ["Ca(2+)", "calcium"],
    "chloride":      ["chloride"],
    "sodium":        ["sodium"],
    "ammonium":      ["ammonium"],
    "manganese":     ["Mn(2+)"],
    "zinc":          ["Zn(2+)"],
    "copper":        ["Cu2(+)", "Cu(2+)"],
    "iron":          ["iron(2+)"],
    "iodide":        ["iodide"],
    "cobalt":        ["cobalt(2+)", "cobalt"],
    "borate":        ["borate", "boric acid"],
    "molybdate":     ["molybdate"],
    "ergosterol":    ["ergosterol"],
    "oleate":        ["oleate", "(9Z)-oleate"],
    "myo_inositol":  ["myo-inositol"],
    "pantothenate":  ["(R)-pantothenate", "pantothenate"],
    "nicotinate":    ["nicotinate"],
    "pyridoxine":    ["pyridoxine"],
    "thiamine":      ["thiamine", "thiamine(1+)"],
    "biotin":        ["biotin", "biotinate"],
    "alanine":       ["L-alanine", "L-alanine zwitterion"],
    "arginine":      ["L-arginine", "L-argininium"],
    "aspartate":     ["L-aspartate"],
    "cysteine":      ["L-cysteine", "L-cysteine zwitterion"],
    "glutamate":     ["L-glutamate"],
    "glutamine":     ["L-glutamine", "L-glutamine zwitterion"],
    "glycine":       ["L-glycine", "glycine", "glycine zwitterion"],
    "histidine":     ["L-histidine", "L-histidine zwitterion"],
    "isoleucine":    ["L-isoleucine", "L-isoleucine zwitterion"],
    "leucine":       ["L-leucine", "L-leucine zwitterion"],
    "lysine":        ["L-lysine", "L-lysinium"],
    "methionine":    ["L-methionine", "L-methionine zwitterion"],
    "phenylalanine": ["L-phenylalanine", "L-phenylalanine zwitterion"],
    "serine":        ["L-serine", "L-serine zwitterion"],
    "threonine":     ["L-threonine", "L-threonine zwitterion"],
    "tryptophan":    ["L-tryptophan", "L-tryptophan zwitterion"],
    "tyrosine":      ["L-tyrosine", "L-tyrosine zwitterion"],
    "valine":        ["L-valine", "L-valine zwitterion"],
    # nao vem da receita: abertos a parte
    "water":         ["H2O", "water"],
    "proton":        ["H+"],
    "co2":           ["carbon dioxide"],
    "ash":           ["ash 1 g"],
    # so existem no yeast-GEM; fazem parte da lista do anaerobicModel.m
    "lanosterol":            ["lanosterol"],
    "zimosterol":            ["zymosterol"],
    "14_demetil_lanosterol": ["14-demethyllanosterol"],
    "palmitoleate":          ["palmitoleate", "(9Z)-hexadecenoate"],
}

# Camada basal: nao sao nutrientes limitantes, ou nao existem na receita.
# O FERRO nao consta da seccao 2.2 do Evers, e o MS300 classico tambem nao o
# tem. E indispensavel nos quatro modelos (T8 do obj1a e do obj1b), por isso
# sem ele nenhum cresce. Fica aberto e declarado, e nao derivado da receita.
BASAIS = ("water", "proton", "co2", "iron")

# Camada de auxotrofia anaerobia e cofactores. Abertos sem tecto, como faz o
# protocolo anaerobio do proprio yeast-GEM (`anaerobicModel.m`, que abre
# ergosterol, lanosterol, zimosterol, 14-demetil-lanosterol, palmitoleato,
# oleato, nicotinato e pantotenato).
#
# Porque nao derivados da receita: a biomassa destes modelos tem fraccoes fixas
# de esterol e de cofactores, e uma celula anaerobia limitada em esterol
# economiza -- desce a fraccao de esterol ate uma decima do normal. Com o tecto
# derivado dos 3 mg/L de ergosterol do MS300 o yeast-GEM fica em mu = 0.007, e
# uma fermentacao com este mosto faz-se na pratica. O tecto estaria a medir a
# rigidez da equacao de biomassa, nao a composicao do mosto.
SUPLEMENTOS = ("ergosterol", "lanosterol", "zimosterol", "14_demetil_lanosterol",
               "palmitoleate", "oleate")

# Tecto dos suplementos, mmol/gDW/h. Medido, nao escolhido: e o patamar em que
# o fluxo de cada esterol/UFA deixa de subir quando se relaxa o tecto
# (ergosterol 0.0078, palmitoleato 0.0196, oleato 0.0065 no yeast-GEM a
# mu = 0.21). Posto acima da procura e abaixo do que os tornaria fonte de
# carbono. As VITAMINAS ficam de fora desta camada: sao derivadas da receita,
# porque abertas sem tecto o modelo come-as como carbono (a piridoxina vai a
# 1.7 mmol/gDW/h com tecto de 1000, 30000x a procura da biomassa).
Q_SUPLEMENTO = 0.02


def _norm(s):
    return " ".join((s or "").strip().lower().split())


def construir_mapa(model):
    """{especie: id da troca | None} para este modelo, pelos SINONIMOS."""
    ex = {}
    for r in model.exchanges:
        ex.setdefault(_norm(next(iter(r.metabolites)).name), r.id)
    return {k: next((ex[_norm(n)] for n in v if _norm(n) in ex), None)
            for k, v in SINONIMOS.items()}


# ═══════════════════════════════════════════════════ 5 · aplicar ao modelo ══
def aplicar_meio_mosto(model, mapa=None, mu=MU_REF, x_final=X_REF,
                       eficiencia=None, q_o2=0.0, ash=None,
                       q_suplemento=None, verbose=True):
    """Aplica o MS300 a um modelo, em camadas. Altera in place: usa `with model:`.

    Camada 1   fecha toda a captacao (secrecao fica livre)
    Camada 2   BASAIS sem tecto: agua, protao, CO2 e ferro (o MS300 nao tem ferro)
    Camada 2b  SUPLEMENTOS sem tecto: esteroides, UFA e vitaminas, como no
               protocolo anaerobio do yeast-GEM. Ver a nota em SUPLEMENTOS
    Camada 3   `ash` para os modelos CarveFungi, que nao tem trocas minerais
               separadas. Aberto e declarado: no eixo mineral estes modelos nao
               sao comparaveis com o yeast-GEM
    Camada 4   bounds derivados da composicao: acucares, acidos organicos,
               fontes de azoto e minerais. E esta a camada que carrega a
               informacao do mosto, e e nela que os bounds prendem
    Camada 5   oxigenio no nivel pedido

    Devolve (aplicados, ausentes). `ausentes` e o que a receita tem e o modelo
    nao: sai no relatorio e nao passa em silencio.
    """
    mapa = construir_mapa(model) if mapa is None else mapa
    bounds = bounds_do_meio(mu=mu, x_final=x_final, eficiencia=eficiencia)
    q_suplemento = Q_SUPLEMENTO if q_suplemento is None else q_suplemento
    # `ash` nao e nutriente, e um lugar na estequiometria: 1 g de minerais por
    # 1/0.048 gDW. Aberto sem tecto deixa de haver limite mineral e o mu
    # dispara (2.60 h-1 nos tres CarveFungi). Amarra-se a biomassa que serve.
    ash = ASH_POR_GDW * mu * 1.25 if ash is None else ash

    for r in model.exchanges:                                   # 1
        r.lower_bound = 0.0

    for k in BASAIS:                                            # 2
        if mapa.get(k):
            model.reactions.get_by_id(mapa[k]).lower_bound = -1000.0

    for k in SUPLEMENTOS:                                       # 2b
        if mapa.get(k):
            model.reactions.get_by_id(mapa[k]).lower_bound = -abs(q_suplemento)

    if mapa.get("ash") and ash:                                 # 3
        model.reactions.get_by_id(mapa["ash"]).lower_bound = -abs(ash)

    aplicados, ausentes = {}, []                                # 4
    for especie, b in bounds.items():
        if especie in BASAIS or especie in SUPLEMENTOS:
            continue                       # ja abertos nas camadas 2 e 2b
        rid = mapa.get(especie)
        if rid is None:
            ausentes.append(especie)
            continue
        model.reactions.get_by_id(rid).lower_bound = b
        aplicados[especie] = (rid, b)

    o2 = next((r.id for r in model.exchanges                    # 5
               if _norm(next(iter(r.metabolites)).name) in ("oxygen", "dioxygen")), None)
    if o2:
        model.reactions.get_by_id(o2).lower_bound = -abs(float(q_o2))

    if verbose:
        print(f"  aplicados: {len(aplicados)} | ausentes do modelo: "
              f"{len(ausentes)} -> {ausentes}")
        print(f"  O2 = {q_o2} mmol/gDW/h" + ("" if o2 else "  (troca de O2 nao encontrada)"))
    return aplicados, ausentes


def tabela_do_meio(model, mapa=None, mu=MU_REF, x_final=X_REF, eficiencia=None):
    """A tabela auditavel: concentracao, bound, e a verificacao de volta."""
    mapa = construir_mapa(model) if mapa is None else mapa
    comp = composicao_mmol_L()
    eficiencia = eficiencia or {}
    linhas = []
    for especie, mmol in sorted(comp.items(), key=lambda x: -x[1]):
        eff = eficiencia.get(especie, EFICIENCIA_OMISSAO)
        b = calcular_bound(mmol, mu, x_final, eff)
        linhas.append({
            "especie": especie,
            "mmol_L": round(mmol, 4),
            "troca": mapa.get(especie) or "(ausente)",
            "eficiencia": eff,
            "bound": round(b, 5),
            "verifica_mmol_L": round(verificar_bound(b, mu, x_final), 4),
            "N_atomos": N_POR_MOLECULA.get(especie, 0),
        })
    return pd.DataFrame(linhas)
