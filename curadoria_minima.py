"""
curadoria_minima.py -- catalogo de curadoria reconstruido de raiz a partir dos
testes T1-T14 dos notebooks `obj1a_gem_sc` e `obj1b_gem_ns` (2026-09-24).

Nao tem relacao com `curadoria.py`: nenhuma edicao foi herdada, todas foram
derivadas de um achado medido nestes dois notebooks.

REGRAS
------
1. Uma edicao so existe se um teste a mostrar. O campo `achado` diz qual, e a
   edicao sem achado nao entra no ficheiro.
2. O campo `arbitro` e a referencia que decide o sentido da correccao. Sem
   arbitro nao ha edicao: uma alteracao sem literatura e uma opiniao.
3. O campo `variante` diz onde se aplica. "anaerobia" nao altera o modelo
   aerobio; "ambas" altera os dois.
4. `aplicar_por_omissao = False` marca as edicoes que criam capacidade nova e
   precisam de uma decisao explicita (verificacao no genoma) antes de entrar.

Cada funcao devolve uma lista de linhas de log com o antes e o depois, para a
seccao 11.1 dos notebooks poder imprimir o que mudou sem reconstruir nada.
"""

from dataclasses import dataclass, field
from typing import Callable, List, Dict, Tuple

import cobra


# ---------------------------------------------------------------------------
# infraestrutura
# ---------------------------------------------------------------------------
@dataclass
class Edicao:
    id: str
    titulo: str
    achado: str          # o teste que a autoriza
    arbitro: str         # a referencia que decide o sentido
    variante: str        # "base" | "anaerobia" | "ambas"
    organismos: Tuple[str, ...]
    accao: Callable
    aplicar_por_omissao: bool = True


CATALOGO: Dict[str, Edicao] = {}


def _reg(e: Edicao):
    CATALOGO[e.id] = e
    return e


def _log(op, rid, antes, depois, nota=""):
    return {"op": op, "reaccao": rid, "antes": antes, "depois": depois, "nota": nota}


def _protao(m, compartimento):
    """O metabolito H+ de um compartimento, pela formula e nao pelo nome."""
    for s in m.metabolites:
        if s.formula == "H" and (s.compartment or "").lower() == compartimento.lower():
            return s
    return None


def _reequilibrar_carga(m, rid):
    """Repoe o balanco de carga de uma reaccao ajustando o H+ do compartimento
    maioritario. E o que o protocolo anaerobio do yeast-GEM faz depois de mexer
    em r_4598: repor o coeficiente do protao pela soma das cargas."""
    r = m.reactions.get_by_id(rid)
    desvio = r.check_mass_balance().get("charge", 0.0)
    if abs(desvio) < 1e-12:
        return None
    comp = max({(s.compartment or "c") for s in r.metabolites},
               key=lambda c: sum(1 for s in r.metabolites if s.compartment == c))
    h = _protao(m, comp)
    if h is None:
        return None
    antes = r.reaction
    r.add_metabolites({h: -desvio})
    return _log("carga", rid, antes, r.reaction, f"desvio de carga {desvio:+.6g} reposto em {h.id}")


def edicoes_de(organismo, variante="ambas", so_por_omissao=True):
    """Lista as edicoes aplicaveis a um organismo ('Sc','Td','Mp','Lt')."""
    out = []
    for e in CATALOGO.values():
        if organismo not in e.organismos:
            continue
        if variante != "ambas" and e.variante not in ("ambas", variante):
            continue
        if so_por_omissao and not e.aplicar_por_omissao:
            continue
        out.append(e)
    return out


def aplicar(m, organismo, ids, verbose=True):
    """Aplica as edicoes indicadas, por id, e devolve o log.
    `ids` e explicito de proposito: nunca se aplica um catalogo inteiro."""
    log = []
    for eid in ids:
        e = CATALOGO[eid]
        if organismo not in e.organismos:
            raise ValueError(f"{eid} nao se aplica a {organismo}")
        linhas = e.accao(m, organismo) or []
        for l in linhas:
            l["edicao"] = eid
        log.extend(linhas)
        if verbose:
            print(f"[{eid}] {e.titulo}  ({len(linhas)} operacoes)")
    return log


# ===========================================================================
# organismos sem variante anaerobia
# ===========================================================================
SEM_ANAEROBIOSE = {
    "Mp": (
        "M. pulcherrima nao se simula em anaerobiose estrita. Decisao de "
        "2026-09-24, registada na seccao 11.2 do obj1b.\n"
        "Razao: e fracamente fermentativa e dependente de O2. No ficheiro por "
        "curar nao crescia a O2 = 0, e as edicoes NS1, NS2b e NS3 davam-lhe "
        "crescimento anaerobio por remocao de artefactos do template. O "
        "resultado era que passava a ser o melhor fermentador anaerobio dos "
        "quatro modelos, o contrario do que a especie faz. O bloqueio era "
        "artefacto, mas a ausencia de crescimento anaerobio nao e: e fenotipo. "
        "Simula-se sempre com O2 declarado (ver Q_O2_DECLARADO)."
    ),
}

# Niveis de O2 declarados para os organismos sem variante anaerobia,
# em mmol/gDW/h. PROVISORIOS: tem de ser amarrados ao protocolo de arejamento
# dos ensaios. Enquanto nao houver protocolo, reporta-se o intervalo e nao um
# ponto (ver a curva de O2 na seccao 11.2 do obj1b).
Q_O2_DECLARADO = {
    "Mp": {
        "microaerobiose": 0.5,   # O2 que entra na inoculacao e nas remontagens
        "arejamento": 2.0,       # arejamento controlado, o regime em que a Mp
                                 # desvia carbono para acetato, CO2 e biomassa
    },
}


# ===========================================================================
# obj1a -- Saccharomyces cerevisiae (yeast-GEM)
# ===========================================================================
def _sc1(m, org):
    """Heme a fora da pseudo-reaccao de cofactores, na variante anaerobia."""
    linhas = []
    r = m.reactions.get_by_id("r_4598")
    heme = m.metabolites.get_by_id("s_3714")
    if heme not in r.metabolites:
        return [_log("nada", "r_4598", "-", "-", "heme a ja nao esta em r_4598")]
    coef = r.metabolites[heme]
    antes = r.reaction
    r.subtract_metabolites({heme: coef})
    linhas.append(_log("remover", "r_4598", antes, r.reaction,
                       f"heme a, coef {coef:+.3g} mmol/gDW "
                       f"(~{abs(coef) * (heme.formula_weight or 0):.2e} mg/gDW)"))
    l = _reequilibrar_carga(m, "r_4598")
    if l:
        linhas.append(l)
    return linhas


_reg(Edicao(
    id="SC1",
    titulo="Heme a fora dos cofactores da biomassa anaerobia",
    achado="T11 + T12: s_3714 em r_4598 (coef 1e-06) e o unico termo que poe "
           "mu = 0 a O2 = 0 com esterois e UFA fornecidos",
    arbitro="yeast-GEM, code/otherChanges/anaerobicModel.m, que poe o coeficiente "
            "de s_3714 em r_4598 a zero e reequilibra o H+; biossintese do heme "
            "em levedura e dependente de O2 em Hem13, Hem14 e Cox15",
    variante="anaerobia",
    organismos=("Sc",),
    accao=_sc1,
))


# ===========================================================================
# obj1b -- nao-Saccharomyces (CarveFungi)
# ===========================================================================
def _ns0(m, org):
    """Formulas do par C18 e a estequiometria de protoes que dependia delas."""
    linhas = []
    alvo = {"ocda": ("C18H35O2", -1), "9Zocda": ("C18H33O2", -1)}
    for base, (formula, carga) in alvo.items():
        for s in [x for x in m.metabolites
                  if x.id == base or x.id.startswith(base + "_")]:
            if s.formula == formula:
                continue
            antes = f"{s.formula} carga={s.charge}"
            s.formula, s.charge = formula, carga
            linhas.append(_log("formula", s.id, antes, f"{formula} carga={carga}",
                               "estearato C18:0 e oleato C18:1 estavam com a mesma formula"))
    # as reaccoes cuja estequiometria de protoes tinha sido escrita contra as
    # formulas erradas: o desvio medido diz quantos protoes faltam ou sobram
    for r in m.reactions:
        if r.boundary:
            continue
        if not any(s.id.startswith(("ocda", "9Zocda", "stcoa", "9Zocdcoa"))
                   for s in r.metabolites):
            continue
        d = {k: v for k, v in r.check_mass_balance().items()
             if k != "charge" and abs(v) > 1e-6}
        if set(d) - {"H"}:
            continue                      # desvio que nao e so de protoes: nao se toca
        dH = d.get("H", 0.0)
        if abs(dH) < 1e-6:
            continue
        comp = (r.id.rsplit("_", 1)[-1] or "c").lower()
        h = _protao(m, comp)
        if h is None:
            continue
        antes = r.reaction
        r.add_metabolites({h: -dH})
        resto = {k: v for k, v in r.check_mass_balance().items()
                 if k != "charge" and abs(v) > 1e-6}
        linhas.append(_log("estequiometria", r.id, antes, r.reaction,
                           f"H {dH:+.0f}; {'equilibrada' if not resto else resto}"))
    return linhas


_reg(Edicao(
    id="NS0",
    titulo="Formulas do par C18 e os protoes que dependiam delas",
    achado="T5.1 (celulas 18-21 do obj1b): ocda (estearato) e 9Zocda (oleato) "
           "estao ambos anotados C18H34O2 nos tres modelos. Com as formulas "
           "corrigidas, seis reaccoes ficam com desvio de 1 protao",
    arbitro="Quimica dos anioes de acido gordo: C18:0 e C18H35O2- e C18:1 e "
            "C18H33O2-; a diferenca de 2 H e o que distingue saturado de "
            "insaturado. Com as duas iguais, a dessaturacao fica elementarmente "
            "nula e o modelo fabrica insaturados sem O2 e sem dessaturase. "
            "Os pares C16 e C14 do mesmo ficheiro estao certos, logo e erro "
            "isolado e nao convencao",
    variante="ambas",
    organismos=("Td", "Mp", "Lt"),
    accao=_ns0,
))


def _ns1(m, org):
    """Reoxidacao anaerobia do quinol: fumarato-redutase explicita.

    Nao se torna UF00664_M reversivel. Adiciona-se uma reaccao nova, so no
    sentido redutor, para a actividade ficar visivel no log e nao se criar um
    par SDH/FRD reversivel a competir com a cadeia respiratoria."""
    if "FRD_ANA_M" in m.reactions:
        return [_log("nada", "FRD_ANA_M", "-", "-", "ja existe")]
    r = cobra.Reaction("FRD_ANA_M")
    r.name = "fumarate reductase (quinol) -- actividade anaerobia"
    r.bounds = (0.0, 1000.0)
    m.add_reactions([r])
    r.add_metabolites({m.metabolites.get_by_id("fum_m"): -1.0,
                       m.metabolites.get_by_id("q6h2_m"): -1.0,
                       m.metabolites.get_by_id("succ_m"): 1.0,
                       m.metabolites.get_by_id("q6_m"): 1.0})
    r.notes["curadoria"] = ("NS1: sem via de reoxidacao anaerobia do quinol, a DHOD "
                            "UF00665_CM nao roda sem O2 e caem as pirimidinas. GPR por "
                            "preencher: confirmar homologo de OSM1/FRD1 no genoma.")
    return [_log("adicionar", "FRD_ANA_M", "(nao existia)", r.reaction,
                 "GPR por preencher -- ver nota da reaccao")]


_reg(Edicao(
    id="NS1",
    titulo="Reoxidacao anaerobia do quinol (fumarato-redutase)",
    achado="T11: dna_1g_c e rna_1g_c so se fabricam com O2. UF00664_M "
           "(succ + q6 --> fum + q6h2) e irreversivel e nao ha outro consumidor "
           "anaerobio de q6h2_m, logo a DHOD UF00665_CM nao roda",
    arbitro="Arikawa et al. 1998 FEMS Microbiol Lett 165:111 e Camarasa et al. "
            "2007 Yeast 24:391 (Osm1/Frd1 sao necessarias para crescimento "
            "anaerobio em S. cerevisiae); Bouwknegt et al. 2021 Fungal Biol "
            "Biotechnol 8:10 (DHOD de classe II suporta anaerobiose)",
    variante="anaerobia",
    organismos=("Td", "Mp", "Lt"),
    accao=_ns1,
    aplicar_por_omissao=False,      # cria capacidade: exige verificacao no genoma
))


def _tirar_da_biomassa(m, ids):
    bio = m.reactions.get_by_id("BIOMASS")
    linhas = []
    for mid in ids:
        if mid not in m.metabolites:
            continue
        s = m.metabolites.get_by_id(mid)
        if s not in bio.metabolites:
            continue
        coef = bio.metabolites[s]
        antes = f"{s.id} com coef {coef:+.4g}"
        bio.subtract_metabolites({s: coef})
        try:
            mg = abs(coef) * (s.formula_weight or 0.0)
        except Exception:
            mg = 0.0
        linhas.append(_log("remover", "BIOMASS", antes, f"{s.id} fora da equacao",
                           f"massa retirada ~{mg:.2e} mg/gDW"))
    return linhas


def _ns2a(m, org):
    return _tirar_da_biomassa(m, ["retn_c"])


_reg(Edicao(
    id="NS2a",
    titulo="Retinoato fora da biomassa",
    achado="T11: retn_c so se fabrica com O2 e e um dos termos que bloqueiam "
           "a anaerobiose. T1/T6: identico nos tres modelos, logo vem do template",
    arbitro="Nao ha via de carotenoides nem de retinoides descrita em "
            "Saccharomycetaceae nem em Metschnikowia; o pigmento da M. pulcherrima "
            "e a pulcherrimina, quelato de ferro derivado da leucina. O yeast-GEM, "
            "curado a mao, nao tem retinoato na biomassa",
    variante="ambas",
    organismos=("Td", "Mp", "Lt"),
    accao=_ns2a,
))


def _ns2b(m, org):
    return _tirar_da_biomassa(m, ["ficytc_c", "focytc_c", "q6_m"])


_reg(Edicao(
    id="NS2b",
    titulo="Cofactores respiratorios fora da biomassa anaerobia",
    achado="T11: ficytc_c, focytc_c e q6_m so se fabricam com O2 (coeficientes "
           "de 1e-04 cada). T12: massa desprezavel, efeito binario em mu",
    arbitro="Citocromo c depende de heme (biossintese dependente de O2 em Hem13/"
            "Hem14); a biossintese da ubiquinona tem hidroxilacoes dependentes de "
            "O2 (Coq6, Coq7). Mesmo criterio do protocolo anaerobio do yeast-GEM "
            "para o heme a",
    variante="anaerobia",
    organismos=("Td", "Mp", "Lt"),
    accao=_ns2b,
))


def _ns3(m, org):
    linhas = []
    r = m.reactions.get_by_id("UF02415_R")
    z = m.metabolites.get_by_id("zymst_1g_r")
    e = m.metabolites.get_by_id("ergst_1g_r")
    if z in r.metabolites:
        antes = r.reaction
        coef_z = r.metabolites[z]
        r.subtract_metabolites({z: coef_z})
        r.add_metabolites({e: coef_z})          # a escala de 1 g mantem-se
        linhas.append(_log("estequiometria", "UF02415_R", antes, r.reaction,
                           f"{abs(coef_z):.2f} de zimosterol passam para ergosterol; "
                           f"soma dos coeficientes mantem-se em 1.00"))
    if "UF02416_R" in m.reactions:
        rc = m.reactions.get_by_id("UF02416_R")
        antes = str(rc.bounds)
        rc.bounds = (0.0, 0.0)
        linhas.append(_log("bounds", "UF02416_R", antes, str(rc.bounds),
                           "rota do colesterol: nao e esterol de levedura"))
    return linhas


_reg(Edicao(
    id="NS3",
    titulo="Pool de esterol a 100 % ergosterol na variante anaerobia",
    achado="T11: sterol_1g_r so se fabrica com O2. UF02415_R exige 0.04 de "
           "zymst_1g_r, que nao tem troca nem sintese anaerobia: 4 % do pool "
           "bloqueiam os 100 %",
    arbitro="Erg1, Erg11 e Erg25 consomem O2 molecular (auxotrofia anaerobia "
            "classica, Andreasen & Stier 1953/1954). O protocolo anaerobio do "
            "yeast-GEM abre a captacao de zimosterol, lanosterol e "
            "14-demetil-lanosterol precisamente por causa destes termos "
            "intermedios; nestes ficheiros nao ha esses metabolitos no "
            "extracelular, por isso corrige-se a composicao do pool",
    variante="anaerobia",
    organismos=("Td", "Mp", "Lt"),
    accao=_ns3,
))


# teto de L-lactato da Lt, mmol/gDW/h, para 20 mmol/gDW/h de hexose.
# 9 g/L de acido lactico sobre ~200 g/L de acucar = 0.09 mol lactato / mol hexose.
LT_TECTO_LACTATO = 0.09 * 20.0


def _ns4(m, org):
    linhas = []
    # (a) D-lactato: sem transportador de efluxo conhecido, nos tres
    if "UF02971_E" in m.reactions:
        r = m.reactions.get_by_id("UF02971_E")
        antes = str(r.bounds)
        r.upper_bound = 0.0
        linhas.append(_log("bounds", "UF02971_E", antes, str(r.bounds),
                           "secrecao de D-lactato fechada"))
    # (b) reducao piruvato -> D-lactato dependente de NAD (so a Td a tem)
    for rid in ("UF00106_c", "UF00106_m"):
        if rid in m.reactions:
            r = m.reactions.get_by_id(rid)
            antes = str(r.bounds)
            r.bounds = (0.0, 0.0)
            linhas.append(_log("bounds", rid, antes, str(r.bounds),
                               "reduccao pyr -> D-lactato: sem suporte em levedura"))
    # (c) L-lactato/NAD
    if "UF00103_C" in m.reactions:
        r = m.reactions.get_by_id("UF00103_C")
        antes = str(r.bounds)
        if org == "Lt":
            linhas.append(_log("manter", "UF00103_C", antes, antes,
                               "L. thermotolerans produz acido L-lactico: reaccao mantida"))
        else:
            r.bounds = (0.0, 1000.0)
            linhas.append(_log("bounds", "UF00103_C", antes, str(r.bounds),
                               "so no sentido oxidativo (lactato -> piruvato)"))
    # (d) teto do L-lactato
    if "UF02967_E" in m.reactions:
        r = m.reactions.get_by_id("UF02967_E")
        antes = str(r.bounds)
        r.upper_bound = LT_TECTO_LACTATO if org == "Lt" else 0.0
        linhas.append(_log("bounds", "UF02967_E", antes, str(r.bounds),
                           f"teto de secrecao de L-lactato ({'dados de vinificacao' if org == 'Lt' else 'nao documentado nesta especie'})"))
    return linhas


_reg(Edicao(
    id="NS4",
    titulo="Lactato: sentido oxidativo e tecto de secrecao",
    achado="T13: com o meio de trabalho o carbono vai todo para lactato "
           "(32-35 mmol/gDW/h) e o etanol fica a zero em todos os niveis de O2. "
           "A celula 48 do obj1b mostra que fechar o lactato nao muda mu: sao "
           "optimos alternativos e a pFBA escolhe por numero de passos",
    arbitro="As Saccharomycetaceae nao tem LDH dependente de NAD no sentido "
            "redutor; o metabolismo do lactato faz-se por Cyb2 (EC 1.1.2.3) e "
            "Dld1, ambas oxidativas, e produzir lactato em S. cerevisiae exige "
            "LDH heterologa. Excepcao: L. thermotolerans (Hranilovic et al. 2021; "
            "Vaquero et al. 2022). T. delbrueckii distingue-se por baixa acidez "
            "volatil, nao por lactato (Benito 2018)",
    variante="ambas",
    organismos=("Td", "Mp", "Lt"),
    accao=_ns4,
))


def _ns5(m, org):
    linhas = []
    for rid, nota in (("UF00058_C", "fru + NADH -> sorbitol"),
                      ("UF00049_C", "fru + NADPH -> manitol")):
        if rid in m.reactions:
            r = m.reactions.get_by_id(rid)
            antes = str(r.bounds)
            r.bounds = (0.0, 0.0)
            linhas.append(_log("bounds", rid, antes, str(r.bounds),
                               f"{nota}: escrita so no sentido redutor"))
    # fechar a secrecao: UF00260_C (manitol/NAD) e reversivel e reabre o mesmo
    # dreno por outro caminho. O tecto e no que sai da celula, que e o que se
    # mede no vinho.
    for r in m.exchanges:
        nome = (next(iter(r.metabolites)).name or "").lower()
        if any(p in nome for p in ("mannitol", "sorbitol", "arabinitol", "xylitol")):
            antes = str(r.bounds)
            r.upper_bound = 0.0
            linhas.append(_log("bounds", r.id, antes, str(r.bounds),
                               f"secrecao de {nome} fechada"))
    return linhas


_reg(Edicao(
    id="NS5",
    titulo="Reducao de hexose a poliol fechada",
    achado="T13/T9: a Td secreta 1.52 de sorbitol e a Lt 1.38 no meio de "
           "trabalho; fechado o sorbitol, o dreno passa para o manitol",
    arbitro="O dreno redox anaerobio da levedura vinica e o glicerol, por "
            "Gpd1/Gpd2 e Gpp1/Gpp2; os mutantes gpd1 gpd2 nao crescem sem "
            "oxigenio (Ansell et al. 1997 EMBO J 16:2179). As desidrogenases "
            "homologas em S. cerevisiae (Sor1/Sor2) actuam no sentido oxidativo. "
            "O manitol no vinho e marcador de alteracao por bacterias lacticas "
            "heterofermentativas a reduzir frutose, nao produto de levedura",
    variante="ambas",
    organismos=("Td", "Mp", "Lt"),
    accao=_ns5,
))


def _ns6(m, org, ngam=0.7):
    if "ATPM" in m.reactions:
        return [_log("nada", "ATPM", "-", "-", "ja existe")]
    r = cobra.Reaction("ATPM")
    r.name = "manutencao nao associada ao crescimento (NGAM)"
    r.bounds = (float(ngam), 1000.0)
    m.add_reactions([r])
    r.add_metabolites({m.metabolites.get_by_id("atp_c"): -1.0,
                       m.metabolites.get_by_id("h2o_c"): -1.0,
                       m.metabolites.get_by_id("adp_c"): 1.0,
                       m.metabolites.get_by_id("pi_c"): 1.0,
                       m.metabolites.get_by_id("h_c"): 1.0})
    r.notes["curadoria"] = ("NS6: valor transferido do yeast-GEM (r_4046, 0.7 "
                            "mmol ATP/gDW/h, Nissen et al. 1997) por nao haver "
                            "dados de manutencao para estas tres especies.")
    return [_log("adicionar", "ATPM", "(nao existia)", r.reaction,
                 f"piso de {ngam} mmol ATP/gDW/h, valor transferido da S. cerevisiae")]


_reg(Edicao(
    id="NS6",
    titulo="Manutencao nao associada ao crescimento (NGAM)",
    achado="T1/inspeccao: nenhuma reaccao com lower_bound > 0 envolvendo ATP, e "
           "nao existe ATPM. So ha GAM (60 no BIOMASS)",
    arbitro="yeast-GEM r_4046 = 0.7 mmol ATP/gDW/h (Nissen et al. 1997). Valor "
            "transferido de outra especie: declarar como tal",
    variante="ambas",
    organismos=("Td", "Mp", "Lt"),
    accao=_ns6,
    aplicar_por_omissao=False,      # valor emprestado: exige decisao explicita
))


# ===========================================================================
# variante anaerobia: funcao, nao segundo ficheiro
# ===========================================================================
def modo_anaerobio(m, organismo, o2_rid=None, ids=None, verbose=False,
                   forcar=False):
    """Aplica em memoria as edicoes de variante "anaerobia" e fecha o O2.

    Porque nao um segundo ficheiro: a anaerobiose e uma CONDICAO de simulacao,
    nao uma correccao do modelo. O heme a, os citocromos, a ubiquinona e o
    zimosterol estao certos no modelo aerobio e so deixam de fazer sentido
    quando nao ha O2 para os sintetizar. E tambem o que o proprio yeast-GEM faz:
    `anaerobicModel.m` e uma funcao aplicada ao modelo consenso, e nao um
    ficheiro `yeast-GEM_anaerobic.xml`.

    A regra de uma so fonte de verdade mantem-se: ha um ficheiro curado por
    organismo, e esta funcao e a unica forma autorizada de o por em anaerobiose.

    `ids=None` aplica as edicoes anaerobias marcadas `aplicar_por_omissao`. As
    que criam capacidade nova (NS1, por exemplo) ficam de fora e tem de ser
    passadas a mao, para a decisao ficar escrita no notebook.
    """
    if organismo in SEM_ANAEROBIOSE and not forcar:
        raise ValueError(SEM_ANAEROBIOSE[organismo] +
                         "\nSe for mesmo essa a intencao, passa forcar=True e "
                         "escreve a justificacao no notebook.")
    if ids is None:
        ids = sorted(e.id for e in CATALOGO.values()
                     if organismo in e.organismos and e.variante == "anaerobia"
                     and e.aplicar_por_omissao)
    ids = list(ids)
    log = aplicar(m, organismo, ids, verbose=verbose)
    if o2_rid and o2_rid in m.reactions:
        m.reactions.get_by_id(o2_rid).lower_bound = 0.0
        l = _log("bounds", o2_rid, "captacao aberta", "lb = 0", "O2 fechado")
        l["edicao"] = "modo_anaerobio"
        log.append(l)
    return log


def modo_microaerobio(m, organismo, o2_rid, regime="microaerobiose", q_o2=None):
    """Abre a captacao de O2 no nivel declarado para este organismo.

    Existe para o numero viver num sitio so. Um organismo que nao se simula em
    anaerobiose tem de ser simulado a um nivel de O2 escrito e justificado, e
    nao ao que cada notebook a jusante decidir usar.

    Os valores de `Q_O2_DECLARADO` sao provisorios: reporta o intervalo da curva
    de O2 enquanto nao houver protocolo de arejamento dos ensaios."""
    if q_o2 is None:
        q_o2 = Q_O2_DECLARADO[organismo][regime]
    m.reactions.get_by_id(o2_rid).lower_bound = -abs(float(q_o2))
    return [_log("bounds", o2_rid, "captacao fechada", f"lb = -{abs(q_o2)}",
                 f"{organismo}: regime '{regime}', valor declarado em "
                 f"Q_O2_DECLARADO (provisorio)")]
