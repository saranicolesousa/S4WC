"""
mapeamento.py -- Identidade dos modelos e mapeamento de trocas.

Fonte unica do CROSS_MAP. Ate 2026-09-09 esta tabela vivia dentro do notebook
001; o notebook de comunidade precisava dela e copia-la era criar exactamente o
tipo de divergencia que a curadoria unificada acabou de eliminar. O 001 deve
passar a importar daqui na proxima revisao.

CINCO variantes de modelo, nao quatro: a M. pulcherrima entra em `base` e em
`C4a_todas`, porque o `curadoria.PENDENTES` pede que a diferenca da C4a seja
reportada e nao assumida.
"""
from __future__ import annotations
import re

# --- identidade -------------------------------------------------------------
CODIGO   = {"Scer_yeastGEM": "Sc", "Td_auto": "Td", "Mp_auto": "Mp",
            "Lt_auto": "Lt", "Mp_C4a": "Mp"}
VARIANTE = {"Scer_yeastGEM": "base", "Td_auto": "base", "Mp_auto": "base",
            "Lt_auto": "base", "Mp_C4a": "C4a_todas"}
STEM     = {"Scer_yeastGEM": "yeast-GEM", "Td_auto": "model_delbrueckii",
            "Mp_auto": "model_pulcherrima", "Lt_auto": "model_thermotolerans",
            "Mp_C4a": "model_pulcherrima"}
ROTULO   = {"Scer_yeastGEM": "Sc", "Td_auto": "Td", "Mp_auto": "Mp",
            "Lt_auto": "Lt", "Mp_C4a": "Mp+C4a"}
BIOMASSA = {"Scer_yeastGEM": "r_2111", "Td_auto": "BIOMASS", "Mp_auto": "BIOMASS",
            "Lt_auto": "BIOMASS", "Mp_C4a": "BIOMASS"}
ORDEM = ["Td_auto", "Mp_auto", "Mp_C4a", "Lt_auto", "Scer_yeastGEM"]
# a Mp entra numa comunidade por UMA variante de cada vez
EXCLUSIVOS = [("Mp_auto", "Mp_C4a")]

# --- trocas -----------------------------------------------------------------
_ENTRADAS_SC = {
 "glucose":"r_1714","fructose":"r_1709","ammonium":"r_1654","oxygen":"r_1992",
 "co2":"r_1672","water":"r_2100","proton":"r_1832","phosphate":"r_2005",
 "sulfate":"r_2060","iron":"r_1861","potassium":"r_2020","sodium":"r_2049",
 "chloride":"r_4593","copper":"r_4594","manganese":"r_4595","zinc":"r_4596",
 "magnesium":"r_4597","calcium":"r_4600","ergosterol":"r_1757","oleate":"r_2189",
 "palmitoleate":"r_1994","biotin":"r_1671","nicotinate":"r_1967","thiamine":"r_2067",
 "pyridoxine":"r_2028","pantothenate":"r_1548","inositol":"r_1947","citrate":"r_1687",
 "malate":"r_1552","alanine":"r_1873","arginine":"r_1879","aspartate":"r_1881",
 "cysteine":"r_1883","glutamate":"r_1889","glutamine":"r_1891","glycine":"r_1810",
 "histidine":"r_1893","isoleucine":"r_1897","leucine":"r_1899","lysine":"r_1900",
 "methionine":"r_1902","phenylalanine":"r_1903","serine":"r_1906","threonine":"r_1911",
 "tryptophan":"r_1912","tyrosine":"r_1913","valine":"r_1914"}
_ENTRADAS_AUTO = {
 "glucose":"UF03288_E","fructose":"UF03275_E","ammonium":"UF03376_E","oxygen":"UF03382_E",
 "co2":"UF03227_E","water":"UF03314_E","proton":"UF03474_E","phosphate":"UF02765_E",
 "sulfate":"UF03456_E","ash":"UF02549_E","iron":"UF03268_E","ergosterol":"UF03007_E",
 "oleate":"UF02973_E","palmitoleate":"UF02972_E","biotin":"UF02997_E","nicotinate":"UF03373_E",
 "thiamine":"UF03477_E","pyridoxine":"UF02775_E","pantothenate":"UF02766_E","inositol":"UF03344_E",
 "citrate":"UF02577_E","malate":"UF03510_E","alanine":"UF03176_E","arginine":"UF02548_E",
 "aspartate":"UF02552_E","cysteine":"UF02585_E","glutamate":"UF03292_E","glutamine":"UF02644_E",
 "glycine":"UF03304_E","histidine":"UF02667_E","isoleucine":"UF03339_E","leucine":"UF03351_E",
 "lysine":"UF02705_E","methionine":"UF02999_E","phenylalanine":"UF02761_E","serine":"UF03454_E",
 "threonine":"UF03479_E","tryptophan":"UF02823_E","tyrosine":"UF02826_E","valine":"UF03496_E"}
_SAIDAS_SC = {"ethanol":"r_1761","glycerol":"r_1808","acetate":"r_1634","succinate":"r_2056",
 "L-lactate":"r_1551","D-lactate":"r_1546","pyruvate":"r_2033","alphaketo":"r_1586",
 "acetaldehyde":"r_1631","urea":"r_2091","sorbitol":"r_1712","propanol":"r_4659",
 "ethyl-pyruvate":"r_4737"}
_SAIDAS_AUTO = {"ethanol":"UF03263_E","glycerol":"UF03302_E","acetate":"UF03160_E",
 "succinate":"UF02804_E","L-lactate":"UF02967_E","D-lactate":"UF02971_E","pyruvate":"UF02776_E",
 "alphaketo":"UF03174_E","acetaldehyde":"UF03155_E","urea":"UF03492_E","mannitol":"UF03370_E",
 "sorbitol":"UF03451_E","gluconate":"UF02642_E","meso-butanediol":"UF03104_E",
 "RR-butanediol":"UF03204_E"}
# suplementos que a curacao anaerobia abre. Fora: r_1993 (palmitato) e r_2055
# (estearato), saturados e ausentes das duas receitas.
_ANAEROBIO_SC = ["r_1757","r_1915","r_2106","r_2134","r_2137","r_4780"]
_ANAEROBIO_AUTO = ["UF03007_E"]
_NAO_NUT_SC = ["r_2100","r_1832","r_1672"]
_NAO_NUT_AUTO = ["UF03314_E","UF03474_E","UF03227_E"]

CROSS_MAP = {}
for _a in ORDEM:
    _sc = (_a == "Scer_yeastGEM")
    CROSS_MAP[_a] = {
        "biomass": BIOMASSA[_a],
        "ngam": "r_4046" if _sc else "UF01847_CE",
        "entradas": dict(_ENTRADAS_SC if _sc else _ENTRADAS_AUTO),
        "saidas": dict(_SAIDAS_SC if _sc else _SAIDAS_AUTO),
        "anaerobio": list(_ANAEROBIO_SC if _sc else _ANAEROBIO_AUTO),
        "nao_nutrientes": list(_NAO_NUT_SC if _sc else _NAO_NUT_AUTO),
    }

# --- classes de nutriente ---------------------------------------------------
INORGANICOS = ("phosphate","sulfate","iron","potassium","sodium","chloride","copper",
               "manganese","zinc","magnesium","calcium","iodide","cobalt","borate","molybdate")
VITAMINAS = ("nicotinate","thiamine","biotin","pyridoxine","pantothenate","inositol")
ALIAS_RECEITA = {"DL_malate": "malate", "myo_inositol": "inositol"}
COLESTEROL = ("cholest_e", "cholestester_e")

FORMULA_ESPERADA = {
    "phosphate":"HO4P","sulfate":"O4S","ammonium":"H4N","glucose":"C6H12O6","fructose":"C6H12O6",
    "ethanol":"C2H6O","glycerol":"C3H8O3","acetate":"C2H3O2","succinate":"C4H4O4",
    "L-lactate":"C3H5O3","D-lactate":"C3H5O3","pyruvate":"C3H3O3","acetaldehyde":"C2H4O",
    "urea":"CH4N2O","oleate":"C18H33O2","palmitoleate":"C16H29O2",
    "citrate":"C6H5O7","malate":"C4H4O5"}

def elementos(formula):
    if not formula:
        return {}
    d = {}
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        if el:
            d[el] = d.get(el, 0) + (int(n) if n else 1)
    return d

def consorcios(especies=("Sc", "Td", "Lt"), com_mp=True):
    """Todos os subconjuntos de tamanho >= 2, com a Mp a entrar por UMA variante
    de cada vez. Devolve lista de tuplos de alias."""
    from itertools import combinations
    base = [a for a in ORDEM if CODIGO[a] in especies and CODIGO[a] != "Mp"]
    grupos = [base]
    if com_mp:
        grupos = [base + ["Mp_auto"], base + ["Mp_C4a"]]
    vistos, out = set(), []
    for g in grupos:
        for n in range(2, len(g) + 1):
            for c in combinations(g, n):
                if c not in vistos:
                    vistos.add(c); out.append(c)
    return out


# ---------------------------------------------------------------------------
# TROCAS PARTILHAVEIS ENTRE ORGANISMOS                   (acrescentado 2026-09-11)
# ---------------------------------------------------------------------------
# O espaco extracelular partilhado da comunidade da, por construcao, permissao de
# CAPTACAO a todas as trocas de todos os modelos. As tres CarveFungi tem ~250
# trocas cada, geradas automaticamente por template, e a maioria nao corresponde
# a transportador nenhum. O resultado medido na rev20260909: 986 trocas entre
# organismos, das quais 953 sao artefacto -- etilenoglicol num sentido e
# glicolaldeido no outro a 20.84 mmol/gDW/h, oxaloacetato e fumarato a 16.9,
# ciclos de carbono fechados a atravessar duas membranas.
#
# Esta lista e o que um organismo pode CAPTAR do pool para alem do que o meio
# fornece. A SECRECAO fica sempre aberta: um produto que sai continua a sair, e
# a acumular no liquido, mesmo que ninguem o possa captar.
#
# Criterio de entrada: transportador de membrana plasmatica documentado em
# levedura, ou difusao passiva estabelecida. Nao entra nada por ser plausivel.
#
#   ethanol, acetaldehyde   difusao passiva pela bicamada. O acetaldeido e a
#                           especie de cross-feeding documentada em fermentacao
#                           vinaria em cultura mista.
#   glycerol                Fps1p (canal, Luyten et al. 1995, EMBO J 14:1360) e
#                           Stl1p (simporte activo, Ferreira et al. 2005, Mol
#                           Biol Cell 16:2068).
#   L-lactate, D-lactate,   Jen1p, simporte monocarboxilato/H+. Casal et al.
#   pyruvate, acetate       1999, J Bacteriol 181:2620: celulas crescidas em
#                           lactato tem "an additional permease which transports
#                           lactate, pyruvate, acetate, and propionate".
#   urea                    Dur3p, permease activa de ureia.
#   propanol                alcool superior, difusao passiva. So o yeast-GEM tem
#                           a troca; fica para quando o gap-filling da via de
#                           Ehrlich cobrir as CarveFungi.
#
# FORA, e a ausencia e deliberada:
#   succinate   e produto maior da fermentacao e a SECRECAO fica aberta, mas nao
#               ha permease de dicarboxilato documentada em S. cerevisiae -- o
#               Jen1p e monocarboxilato (Casal 1999, acima) e a razao de a
#               engenharia metabolica importar o Mae1p de S. pombe e precisamente
#               essa. Captar succinato do pool seria uma interaccao inventada.
#   sorbitol, mannitol, butanodiois, gluconato, oxaloacetato, fumarato,
#   intermediarios fosforilados, nucleotidos, tioesteres
#               sem transportador de membrana plasmatica conhecido. Sao o grosso
#               dos 953 artefactos.
#
# Os nutrientes do mosto NAO precisam de estar aqui: o `montar_consorcio_livre`
# junta automaticamente tudo o que a receita fornece (acucares, azoto amoniacal e
# aminoacidico, vitaminas, inorganicos, esterol e UFA, agua, protao, CO2). E por
# isso que a complementaridade lipidica -- que e a interaccao que este projecto
# procura -- continua visivel.
# ---------------------------------------------------------------------------
PARTILHAVEIS_SAIDAS = ("ethanol", "glycerol", "acetate", "L-lactate", "D-lactate",
                       "pyruvate", "acetaldehyde", "urea", "propanol")

# ---------------------------------------------------------------------------
# SEGUNDO CRITERIO: REDOX ANAEROBIO                     (acrescentado 2026-09-11)
#
# Ter transportador nao e o mesmo que poder assimilar. A lista acima responde
# "o composto atravessa a membrana?"; falta responder "e o organismo consegue
# metaboliza-lo sem oxigenio?".
#
# Porque foi preciso. Medido no par Td+Sc a O2 = 0: a Sc captava 12.2
# mmol/gDW/h de etanol do pool e 8.3 de malato, e fechava o redox despejando
# isobutanol (6.16) e 2,3-butanodiol (6.04), 49 mmol C/h em dois compostos que
# nenhum ensaio de vinho mostra a esses niveis. O etanol da comunidade descia
# porque um parceiro o comia, o que nao e uma via para vinho de baixo alcool --
# e um artefacto de reconstrucao.
#
# O criterio, e aplica-se so a CAPTACAO: o catabolismo do composto tem de ser
# redox-neutro ou consumidor de NADH em anaerobiose. Se gerar NADH, so a cadeia
# respiratoria o reoxida, e a O2 = 0 isso nao existe.
#
#   FICAM
#     acetaldehyde  reduzido a etanol pela ADH, CONSOME NADH. E o sumidouro de
#                   equivalentes redutores classico e a especie de cross-feeding
#                   documentada em fermentacao vinaria em cultura mista.
#     pyruvate      PDC -> acetaldeido + CO2 (sem redox) -> etanol (consome
#                   NADH). Redox-neutro a consumidor.
#     urea          fonte de AZOTO, nao de carbono; nao entra no balanco redox.
#
#   SAEM
#     ethanol       ADH2 -> acetaldeido -> ALD -> acetato -> ACS -> acetil-CoA.
#                   Os dois primeiros passos geram NADH e o esqueleto so segue
#                   pelo ciclo do glioxilato e gluconeogenese. Sem respiracao
#                   nao ha como reoxidar. O etanol e o produto final da
#                   fermentacao, nao um substrato dela.
#     glycerol      "glycerol is catabolized in a fully respiratory manner by
#                   ... Saccharomyces cerevisiae" e "no fermentation products
#                   have ever been reported in cultivations using S. cerevisiae
#                   wild-type strains growing on glycerol" (Klein et al. 2019,
#                   Biotechnol Biofuels 12:257). A via nativa L-G3P canaliza
#                   electroes por FADH2 para a cadeia respiratoria.
#     acetate       assimilacao por ACS -> acetil-CoA exige glioxilato e TCA,
#                   ambos dependentes de respiracao para fechar o redox.
#     L/D-lactate   a oxidacao a piruvato faz-se pela Cyb2p (citocromo b2) e
#                   pela Dld1p, ambas ligadas a cadeia respiratoria.
#     propanol      alcool superior, sem via de assimilacao. Fica de fora por
#                   omissao e nao por decisao: o fluxo era 1e-6.
#
# A SECRECAO destes compostos continua toda aberta -- o etanol, o glicerol, o
# acetato e o lactato continuam a sair e a acumular no liquido. O que se fecha e
# o caminho de volta para dentro de outro organismo.
#
# Consequencia que tem de ir para os metodos, e e um resultado e nao uma
# limitacao: em anaerobiose estrita, entre leveduras, o cross-feeding reduz-se a
# azoto, lipidos, vitaminas e a ponte do acetaldeido. Nao ha alimentacao cruzada
# de carbono pelos produtos de fermentacao, porque nenhum deles e assimilavel
# sem oxigenio.
# ---------------------------------------------------------------------------
PARTILHAVEIS_ANAEROBIO = ("acetaldehyde", "pyruvate", "urea")

MODOS_PARTILHA = {"anaerobio": PARTILHAVEIS_ANAEROBIO,   # adoptado
                  "transporte": PARTILHAVEIS_SAIDAS}     # so o criterio de transporte


def rids_partilhaveis(alias, modo="anaerobio"):
    """ids de troca do modelo `alias` que podem ser CAPTADOS do pool.

    modo  "anaerobio"  transportador documentado E assimilavel sem O2 (adoptado)
          "transporte" so o criterio de transportador, para sensibilidade
    """
    if modo not in MODOS_PARTILHA:
        raise ValueError(f"modo {modo!r} nao esta em {sorted(MODOS_PARTILHA)}")
    saidas = CROSS_MAP[alias]["saidas"]
    return [saidas[k] for k in MODOS_PARTILHA[modo] if k in saidas]


def consorcios_c4a(especies=("Sc", "Td", "Lt")):
    """Consórcios com a Mp a entrar SO pela variante C4a.

    Decisao de 2026-09-11: a Mp base e infactivel a O2 = 0 (nao constroi membrana
    em anaerobiose estrita), logo nao e parceiro. Fecha a pendencia da C4a e os
    consorcios passam de 18 para 11."""
    from itertools import combinations
    g = [a for a in ORDEM if CODIGO[a] in especies and CODIGO[a] != "Mp"] + ["Mp_C4a"]
    return [c for n in range(2, len(g) + 1) for c in combinations(g, n)]
