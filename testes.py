"""
testes.py -- testes gerais de um GEM, antes de decidir qualquer curadoria.

Nove testes. Cada um e uma funcao curta que faz uma pergunta e devolve a
resposta. Nenhum altera o modelo: todos correm dentro de `with model:` ou
sobre copias, por isso podem ser corridos em qualquer ordem e quantas vezes
for preciso.

A regra destes notebooks: primeiro mede-se, depois decide-se o que curar.
Nenhuma edicao entra sem que um destes testes mostre o que esta partido.
"""

import pandas as pd
from cobra.medium import minimal_medium
from cobra.flux_analysis import pfba


# ---------------------------------------------------------------------------
# auxiliares
# ---------------------------------------------------------------------------
def _num(v):
    """Converte o resultado de slim_optimize em numero (nan -> 0)."""
    return 0.0 if v is None or v != v else float(v)


def _externo(met):
    return (met.compartment or "").lower() in ("e", "ex", "extracellular")


def abrir(m, rids, q=10.0):
    """Abre as trocas da lista a captacao de q mmol/gDW/h. Altera in place."""
    abertas = []
    for rid in rids:
        if rid in m.reactions:
            m.reactions.get_by_id(rid).lower_bound = -abs(q)
            abertas.append(rid)
    return abertas


def fechar_captacao(m):
    """Fecha a captacao de todas as trocas, deixando a secrecao livre."""
    for r in m.exchanges:
        r.lower_bound = 0.0


# ---------------------------------------------------------------------------
# T1 -- o que ha no ficheiro
# ---------------------------------------------------------------------------
def t1_contar(m):
    comps = sorted({(s.compartment or "?") for s in m.metabolites})
    sem_formula = sum(1 for s in m.metabolites if not s.formula)
    return {
        "modelo": m.id,
        "reaccoes": len(m.reactions),
        "metabolitos": len(m.metabolites),
        "genes": len(m.genes),
        "compartimentos": ", ".join(comps),
        "trocas": len(m.exchanges),
        "metabolitos_sem_formula": sem_formula,
        "reaccoes_sem_GPR": sum(1 for r in m.reactions if not r.gene_reaction_rule),
    }


# ---------------------------------------------------------------------------
# T2 / T3 -- cresce como vem? e com tudo aberto?
# ---------------------------------------------------------------------------
def t2_crescimento(m):
    """mu no meio que vem dentro do ficheiro, e as trocas que ele deixa abertas."""
    abertas = [r for r in m.exchanges if r.lower_bound < 0]
    return {
        "mu_meio_do_ficheiro": round(_num(m.slim_optimize()), 6),
        "trocas_abertas_a_captacao": len(abertas),
        "ids": sorted(r.id for r in abertas),
    }


def t3_tudo_aberto(m, q=1000.0):
    """mu com todas as trocas abertas. Separa 'modelo partido' de 'meio apertado'."""
    with m:
        for r in m.exchanges:
            r.lower_bound = -abs(q)
        return round(_num(m.slim_optimize()), 6)


# ---------------------------------------------------------------------------
# T4 -- sistema fechado: fabrica materia ou energia do nada?
# ---------------------------------------------------------------------------
def _fechar_tudo(m):
    for r in m.boundary:
        r.bounds = (0.0, 0.0)
    for r in m.reactions:
        if r.lower_bound > 0:          # NGAM fixa tornaria o sistema inviavel
            r.lower_bound = 0.0


def t4_sistema_fechado(m, ngam_rid=None):
    """Com tudo fechado, o modelo deve dar mu = 0 e ATP = 0. Qualquer valor
    positivo significa que ha uma fuga, e que TODOS os mu deste modelo estao
    inflacionados (Fritzemeier et al. 2017)."""
    out = {}
    with m:
        _fechar_tudo(m)
        out["biomassa_do_nada"] = round(_num(m.slim_optimize()), 8) + 0.0
    if ngam_rid and ngam_rid in m.reactions:
        with m:
            _fechar_tudo(m)
            m.objective = ngam_rid
            m.objective.direction = "max"
            out["ATP_do_nada"] = round(_num(m.slim_optimize()), 6) + 0.0
    return out


# ---------------------------------------------------------------------------
# T5 -- balancos de massa
# ---------------------------------------------------------------------------
def t5_balancos(m):
    """Reaccoes que nao fecham em elementos. Ignora-se o desvio de carga: o que
    permite fabricar materia e o desequilibrio elementar. As que atravessam a
    membrana plasmatica sao as graves, porque fabricam materia a custa do meio.

    A coluna `sem_formula` separa o desequilibrio real do que nao se consegue
    avaliar: uma reaccao com um metabolito sem formula aparece desequilibrada
    por falta de informacao, e nao por erro de estequiometria."""
    linhas = []
    for r in m.reactions:
        if r.boundary:
            continue
        try:
            d = r.check_mass_balance()
        except Exception:
            continue
        d = {k: round(v, 4) for k, v in d.items() if k != "charge" and abs(v) > 1e-6}
        if not d:
            continue
        comps = sorted({(s.compartment or "?") for s in r.metabolites})
        linhas.append({
            "rxn": r.id,
            "nome": (r.name or "")[:44],
            "compartimentos": ",".join(comps),
            "toca_extracelular": any(_externo(s) for s in r.metabolites),
            "sem_formula": any(not s.formula for s in r.metabolites),
            "desvio": d,
        })
    return pd.DataFrame(linhas)


# ---------------------------------------------------------------------------
# T6 -- a biomassa
# ---------------------------------------------------------------------------
def _pseudo(r):
    """Pseudo-reaccao: sem genes e com varios reagentes (agrega componentes)."""
    return (not r.gene_reaction_rule) and len(r.reactants) > 2


_ENERGIA = ("atp", "h2o", "water", "adp", "phosphate", "h+", "proton")


def _e_energia(met):
    n = (met.name or "").strip().lower()
    return any(n == e or n.startswith(e + " ") for e in _ENERGIA)


def t6_biomassa(m, bio_rid, maxprof=3):
    """Tres perguntas: uma unidade de fluxo vale 1 g? quanto ATP custa crescer?
    que componentes o modelo nao sabe fabricar?

    Sobre a escala, ha duas convencoes e o resumo devolve as duas, porque cada
    modelo usa uma:

      massa_por_formula   soma de |coef| x MW dos reagentes que TEM formula,
                          excluindo o par ATP/H2O da manutencao, que e custo
                          energetico e nao materia. Vale para modelos cuja
                          biomassa lista metabolitos reais.
      soma_coef_1g        soma dos coeficientes dos pseudo-componentes `*_1g_*`
                          (1 g de proteina, 1 g de lipido, ...). Nestes modelos
                          a escala e imposta por construcao e o que tem de dar
                          1 e esta soma, nao a massa por formula.

    Devolve (resumo, tabela_de_componentes). A tabela desce as pseudo-reaccoes
    ate aos componentes reais, porque a biomassa e partida em protein, lipid,
    carbohydrate, RNA, DNA, cofactor e ion."""
    bio = m.reactions.get_by_id(bio_rid)

    # --- escala -------------------------------------------------------------
    massa, opacos, coef_opacos, soma_1g = 0.0, 0, 0.0, 0.0
    for met, coef in bio.metabolites.items():
        if coef >= 0:
            continue
        if "_1g" in met.id.lower():
            soma_1g += abs(coef)
            continue
        if _e_energia(met):
            continue
        try:
            mw = met.formula_weight
        except Exception:
            mw = None
        if not mw or "R" in (met.formula or ""):
            opacos += 1
            coef_opacos += abs(coef)
            continue
        massa += abs(coef) * mw / 1000.0

    # --- GAM: ATP consumido por unidade de biomassa -------------------------
    gam = None
    for met, coef in bio.metabolites.items():
        if coef < 0 and (met.name or "").strip().lower().startswith("atp"):
            gam = abs(coef)
            break

    resumo = {
        "reaccao": bio_rid,
        "n_reagentes": len(bio.reactants),
        "massa_por_formula_g": round(massa, 4),
        "componentes_sem_massa": opacos,
        "soma_coef_componentes_sem_massa": round(coef_opacos, 4),
        "soma_coef_1g": round(soma_1g, 4),
        "GAM_mmol_ATP_por_gDW": gam,
    }
    if soma_1g > 0:
        resumo["escala"] = "pseudo-componentes _1g"
        resumo["desvio_de_1g_pct"] = round((soma_1g - 1.0) * 100, 1)
    elif opacos:
        resumo["escala"] = ("imposta por construcao: %d pseudo-componentes sem massa "
                            "calculavel. Verificar um nivel abaixo." % opacos)
        resumo["desvio_de_1g_pct"] = None
    else:
        resumo["escala"] = "massa por formula"
        resumo["desvio_de_1g_pct"] = round((massa - 1.0) * 100, 1)

    # --- folhas: desce as pseudo-reaccoes -----------------------------------
    folhas, fila, vistos = [], [(s, bio_rid, 0) for s in bio.reactants], set()
    while fila:
        met, origem, prof = fila.pop(0)
        if met.id in vistos:
            continue
        vistos.add(met.id)
        produtoras = [r for r in met.reactions if r.metabolites[met] > 0 and not r.boundary]
        if prof < maxprof and len(produtoras) == 1 and _pseudo(produtoras[0]):
            for s in produtoras[0].reactants:
                fila.append((s, produtoras[0].id, prof + 1))
            continue
        sinteses = [r for r in produtoras
                    if len({(s.compartment or "?") for s in r.metabolites}) == 1
                    or len(r.metabolites) > 2]
        folhas.append({
            "metabolito": met.id,
            "nome": (met.name or "")[:36],
            "formula": met.formula,
            "vem_de": origem,
            "n_produtoras": len(produtoras),
            "n_sinteses": len(sinteses),
            "so_por_importacao": len(sinteses) == 0,
            "tem_troca": any(r.boundary for r in met.reactions),
        })
    return resumo, pd.DataFrame(folhas)


# ---------------------------------------------------------------------------
# T7 -- meio minimo: de que e que o modelo precisa mesmo?
# ---------------------------------------------------------------------------
def t7_meio_minimo(m, fraccao=0.1, abrir_tudo=True, minimizar_componentes=False,
                   fontes_de_C=None, q_fonte=10.0, manter_abertas=None):
    """Menor conjunto de trocas que sustenta `fraccao` do mu maximo.

    `abrir_tudo=True` parte de todas as trocas abertas, para responder a "de que
    e que este organismo precisa" e nao a "o que e que este ficheiro deixa
    entrar". `minimizar_componentes=True` minimiza o NUMERO de compostos (MILP,
    lento); por omissao minimiza o fluxo total de importacao (LP, rapido).

    `fontes_de_C` e o argumento que torna a resposta util para este projecto.
    Sem ele, o LP escolhe a fonte de carbono mais densa que encontrar (tipicamente
    um oligossacarido), porque assim minimiza o fluxo total -- resposta correcta
    para a pergunta errada. Com ele, fecha-se a captacao de TODOS os compostos com
    carbono excepto os desta lista, limitados a `q_fonte`, e a pergunta passa a ser
    "com este acucar, de que mais e que o modelo precisa"."""
    with m:
        if abrir_tudo:
            for r in m.exchanges:
                r.lower_bound = -1000.0
        fechadas_sem_formula = []
        if fontes_de_C:
            manter = set(fontes_de_C) | set(manter_abertas or ())
            for r in m.exchanges:
                if r.id in manter:
                    continue
                met = next(iter(r.metabolites))
                if not met.formula:
                    # sem formula nao se consegue excluir que traga carbono;
                    # fecha-se e regista-se, para a decisao ser visivel
                    r.lower_bound = 0.0
                    fechadas_sem_formula.append(f"{r.id} ({met.name})")
                elif "C" in (met.elements or {}):
                    r.lower_bound = 0.0
            for rid in fontes_de_C:
                if rid in m.reactions:
                    m.reactions.get_by_id(rid).lower_bound = -abs(q_fonte)
        if fechadas_sem_formula:
            print(f"   {len(fechadas_sem_formula)} trocas fechadas por o metabolito "
                  f"nao ter formula: {', '.join(fechadas_sem_formula[:6])}"
                  f"{' ...' if len(fechadas_sem_formula) > 6 else ''}")
        mu = _num(m.slim_optimize())
        if mu < 1e-9:
            return None, mu
        meio = minimal_medium(m, mu * fraccao,
                              minimize_components=minimizar_componentes)
    return meio, mu


def descrever_meio(m, meio):
    """Da nome e formula a cada troca do meio minimo."""
    if meio is None:
        return pd.DataFrame()
    linhas = []
    for rid, fluxo in meio.items():
        r = m.reactions.get_by_id(rid)
        met = next(iter(r.metabolites))
        linhas.append({"troca": rid, "metabolito": (met.name or "")[:40],
                       "formula": met.formula, "captacao": round(float(fluxo), 4)})
    return pd.DataFrame(linhas).sort_values("captacao", ascending=False)


# ---------------------------------------------------------------------------
# T8 -- que compostos sao indispensaveis
# ---------------------------------------------------------------------------
def t8_essenciais(m, limiar=1e-6):
    """Fecha uma troca do meio de cada vez e mede o que sobra do crescimento.
    Corre sobre o meio que o modelo tem NESTE momento, por isso a resposta
    depende do meio: e uma propriedade do par modelo+meio, nao do modelo."""
    base = _num(m.slim_optimize())
    linhas = []
    for r in m.exchanges:
        if r.lower_bound >= 0:
            continue
        met = next(iter(r.metabolites))
        with m:
            r.lower_bound = 0.0
            mu = _num(m.slim_optimize())
        linhas.append({
            "troca": r.id,
            "metabolito": (met.name or "")[:36],
            "mu_sem_ele": round(mu, 6),
            "perda_pct": round((base - mu) / base * 100, 1) if base > limiar else None,
            "essencial": mu < limiar,
        })
    df = pd.DataFrame(linhas)
    return df.sort_values(["essencial", "perda_pct"], ascending=[False, False]) \
        if len(df) else df


# ---------------------------------------------------------------------------
# T9 -- aerobiose e anaerobiose
# ---------------------------------------------------------------------------
def _objectivo_rid(m):
    from cobra.util.solver import linear_reaction_coefficients
    alvos = list(linear_reaction_coefficients(m))
    return alvos[0].id if alvos else None


def t9_oxigenio(m, o2_rid, niveis=(0.0, 1.0, 5.0, 1000.0), produtos=None,
                parcimonioso=True, bio_rid=None):
    """mu e secrecoes em funcao do oxigenio disponivel, no meio que o modelo
    tem neste momento. O nivel 0 e a pergunta que decide se o modelo serve para
    simular fermentacao: cresce em anaerobiose estrita?

    `parcimonioso=True` usa pFBA: sem isso os valores de secrecao sao um entre
    muitos optimos alternativos e mudam entre execucoes sem que nada tenha
    mudado no modelo."""
    bio_rid = bio_rid or _objectivo_rid(m)
    linhas = []
    for q in niveis:
        with m:
            if o2_rid in m.reactions:
                m.reactions.get_by_id(o2_rid).lower_bound = -abs(float(q))
            if parcimonioso:
                try:
                    sol = pfba(m)
                except Exception:
                    sol = m.optimize()
            else:
                sol = m.optimize()
            ok = sol.status == "optimal"
            # com pFBA o objective_value e o fluxo total, nao o crescimento:
            # o mu le-se sempre no fluxo da reaccao objectivo
            mu = float(sol.fluxes.get(bio_rid, 0.0)) if ok else 0.0
            linha = {"O2_disponivel": q, "mu": round(mu, 6),
                     "cresce": bool(mu > 1e-6), "estado": sol.status}
            for nome, rid in (produtos or {}).items():
                linha[nome] = round(float(sol.fluxes.get(rid, 0.0)), 4) if ok else 0.0
        linhas.append(linha)
    return pd.DataFrame(linhas)

def localizar(m, palavras):
    """Todas as reaccoes que tocam metabolitos cujo nome contem uma das palavras."""
    linhas, vistas = [], set()
    for s in m.metabolites:
        nome = (s.name or "").lower()
        if not any(p in nome for p in palavras):
            continue
        for r in s.reactions:
            if r.id in vistas or r.boundary:
                continue
            vistas.add(r.id)
            linhas.append({"rxn": r.id, "bounds": str(r.bounds),
                           "equacao": r.reaction[:74],
                           "gpr": r.gene_reaction_rule or "(SEM GENE)"})
    return pd.DataFrame(linhas)


# ===========================================================================
# T10 -- ciclos geradores de energia, por transportador de energia
# ===========================================================================
def t10_egc(modelo, pares):
    """Fritzemeier et al. 2017, versao completa. Com TODAS as fronteiras
    fechadas e nenhuma reaccao forcada (lb > 0 ou ub < 0 anulados), adiciona-se
    uma reaccao de dissipacao por transportador de energia e maximiza-se.
    Qualquer valor > 1e-6 e um ciclo gerador de energia.

    O T4 sozinho nao chega: uma fuga de redox nao aparece no teste do ATP se
    faltar um precursor, e um ciclo de protoes nao aparece em nenhum dos dois.

    Usa-se uma COPIA por teste. Adicionar e remover reaccoes dentro de um
    `with modelo:` corrompe o container de variaveis do solver."""
    import cobra
    out = {}
    for nome, (subs, prods) in pares.items():
        m = modelo.copy()
        for r in m.boundary:
            r.bounds = (0.0, 0.0)
        for r in m.reactions:
            if r.lower_bound > 0:
                r.lower_bound = 0.0
            if r.upper_bound < 0:
                r.upper_bound = 0.0
        try:
            mets = {m.metabolites.get_by_id(k): v
                    for k, v in {**subs, **prods}.items()}
        except KeyError:
            out[nome] = None
            continue
        rx = cobra.Reaction("EGC_TESTE")
        rx.bounds = (0.0, 1000.0)
        m.add_reactions([rx])
        rx.add_metabolites(mets)
        m.objective = rx
        v = m.slim_optimize()
        out[nome] = 0.0 if v != v else round(v, 4)
    return out


# pares de dissipacao, por familia de modelo
PARES_SC = {
    "ATP_c":   ({"s_0434": -1, "s_0803": -1}, {"s_0394": 1, "s_1322": 1, "s_0794": 1}),
    "ATP_m":   ({"s_0437": -1, "s_0807": -1}, {"s_0397": 1, "s_1325": 1, "s_0799": 1}),
    "NADH_c":  ({"s_1203": -1}, {"s_1198": 1, "s_0794": 1}),
    "NADPH_c": ({"s_1212": -1}, {"s_1207": 1, "s_0794": 1}),
}

PARES_NS = {
    "ATP_c":    ({"atp_c": -1, "h2o_c": -1}, {"adp_c": 1, "pi_c": 1, "h_c": 1}),
    "NADH_c":   ({"nadh_c": -1}, {"nad_c": 1, "h_c": 1}),
    "NADPH_c":  ({"nadph_c": -1}, {"nadp_c": 1, "h_c": 1}),
    "Q6H2_m":   ({"q6h2_m": -1}, {"q6_m": 1, "h_m": 2}),
    "protao_m": ({"hi_m": -1}, {"h_m": 1}),
    "ACCOA_c":  ({"accoa_c": -1, "h2o_c": -1}, {"ac_c": 1, "coa_c": 1, "h_c": 1}),
}


# ===========================================================================
# T11 -- o que bloqueia a anaerobiose, componente a componente
# ===========================================================================
def t11_bloqueio_anaerobio(m, bio_rid, aplicar_meio, maxprof=2):
    """Para cada reagente da biomassa -- e, recursivamente, para os reagentes
    das pseudo-reaccoes que produzem um componente bloqueado -- a producao
    maxima com e sem O2, no meio de trabalho.

    Sao os `so_com_O2 == True` que explicam mu = 0 em anaerobiose. A descida
    recursiva e o que distingue este teste de olhar so para a biomassa: no
    yeast-GEM o bloqueio esta dois niveis abaixo (r_4041 -> r_4598 -> heme a).

    `aplicar_meio(m, q_o2=...)` e a funcao de meio do proprio notebook."""
    bio = m.reactions.get_by_id(bio_rid)
    linhas, fila, vistos = [], [(s, bio_rid, 0) for s in bio.reactants], set()
    while fila:
        met, origem, prof = fila.pop(0)
        if met.id in vistos:
            continue
        vistos.add(met.id)
        if (met.formula or "") == "O2":
            continue                      # o proprio oxigenio nao e um achado
        v = {}
        for q in (0.0, 2.0):
            with m:
                aplicar_meio(m, q_o2=q)
                m.objective = m.add_boundary(met, type="demand")
                x = m.slim_optimize()
            v[q] = 0.0 if x != x else round(x, 6)
        bloq = (v[0.0] < 1e-9) and (v[2.0] > 1e-9)
        linhas.append({"metabolito": met.id, "nome": (met.name or "")[:30],
                       "vem_de": origem, "nivel": prof,
                       "max_sem_O2": v[0.0], "max_com_O2": v[2.0],
                       "so_com_O2": bloq})
        if bloq and prof < maxprof:
            prods = [r for r in met.reactions
                     if r.metabolites[met] > 0 and not r.boundary
                     and not r.gene_reaction_rule and len(r.reactants) > 2]
            for r in prods:
                for s in r.reactants:
                    fila.append((s, r.id, prof + 1))
    return pd.DataFrame(linhas)


# ===========================================================================
# T12 -- sensibilidade de mu a cada termo residual da biomassa
# ===========================================================================
def t12_termos_residuais(m, rids_pseudo, aplicar_meio, limiar_coef=1e-4):
    """Retira, um de cada vez, cada reagente com coeficiente abaixo do limiar,
    e mede mu com e sem ele, a O2 = 0 e a O2 = 2.

    O que se procura: um termo cuja massa e desprezavel (coluna
    `massa_mg_por_gDW`) mas que muda mu de zero para positivo. Um termo assim
    nao esta a descrever composicao celular, esta a impor um fenotipo -- e, com
    coeficientes da ordem de 1e-6, a impo-lo a distancia da tolerancia do
    solver, o que significa que a resposta pode mudar so por trocar de solver."""
    linhas = []
    for rid in rids_pseudo:
        if rid not in m.reactions:
            continue
        r0 = m.reactions.get_by_id(rid)
        for met, coef in list(r0.metabolites.items()):
            if coef >= 0 or abs(coef) > limiar_coef:
                continue
            res = {}
            for q in (0.0, 2.0):
                with m:
                    aplicar_meio(m, q_o2=q)
                    y = m.slim_optimize()
                    res[f"com_{q}"] = 0.0 if y != y else round(y, 6)
                with m:
                    aplicar_meio(m, q_o2=q)
                    m.reactions.get_by_id(rid).subtract_metabolites({met: coef})
                    x = m.slim_optimize()
                    res[f"sem_{q}"] = 0.0 if x != x else round(x, 6)
            try:
                mw = met.formula_weight or 0.0
            except Exception:
                mw = 0.0
            linhas.append({
                "pseudo": rid, "metabolito": met.id, "nome": (met.name or "")[:28],
                "coef": coef, "massa_mg_por_gDW": round(abs(coef) * mw, 6),
                "mu_O2_0_com": res["com_0.0"], "mu_O2_0_sem": res["sem_0.0"],
                "mu_O2_2_com": res["com_2.0"], "mu_O2_2_sem": res["sem_2.0"],
                "decide_anaerobiose": res["com_0.0"] < 1e-6 <= res["sem_0.0"],
            })
    return pd.DataFrame(linhas)


# ===========================================================================
# T13 -- espectro de produtos por FVA (substitui as colunas de pFBA do T9)
# ===========================================================================
def t13_produtos(m, o2_rid, saidas, niveis=(0.0, 0.5, 2.0, 20.0), fraccao=0.999):
    """Intervalo (minimo, maximo) de cada produto no optimo, por nivel de O2.

    Como ler:
      minimo > 0                -> o produto e obrigatorio naquele optimo
      minimo = 0 e maximo > 0   -> possivel mas alternativo: nao se pode
                                   afirmar nem que o modelo o produz nem que
                                   nao o produz
      maximo = 0                -> o modelo nao consegue produzi-lo

    O T9 com pFBA da um ponto entre muitos optimos, escolhido por numero de
    passos e nao por bioquimica. Nao serve para afirmar o que o modelo
    'produz'."""
    from cobra.flux_analysis import flux_variability_analysis as fva
    rids = [r for r in saidas.values() if r in m.reactions]
    nomes = {v: k for k, v in saidas.items()}
    out = {}
    for q in niveis:
        with m:
            if o2_rid in m.reactions:
                m.reactions.get_by_id(o2_rid).lower_bound = -abs(float(q))
            mu = m.slim_optimize()
            if mu != mu or mu < 1e-9:
                out[f"O2={q:g}"] = {nomes[r]: "sem crescimento" for r in rids}
                continue
            f = fva(m, reaction_list=rids, fraction_of_optimum=fraccao)
        out[f"O2={q:g}"] = {nomes[r]: (round(f.loc[r, "minimum"], 2),
                                       round(f.loc[r, "maximum"], 2)) for r in rids}
    return pd.DataFrame(out)


# ===========================================================================
# T14 -- balanco de carbono e rendimento em etanol
# ===========================================================================
def _carbonos(m, rid):
    met = next(iter(m.reactions.get_by_id(rid).metabolites))
    return (met.elements or {}).get("C", 0)


def t14_carbono(m, bio_rid, etanol_rid=None, C_por_gDW=41.0):
    """Carbono que entra pelas trocas vs carbono que sai pelas trocas mais o
    carbono que foi para biomassa.

    C_por_gDW = 41 mmol C/gDW corresponde a composicao elementar tipica de
    levedura (CH1.8O0.5N0.2, ~25.8 g/C-mol).

    fecho_pct entre 95 e 105 e o aceitavel. Abaixo, ha carbono a desaparecer
    numa reaccao desequilibrada; acima, esta a ser criado. Um pseudo-metabolito
    sem formula (por exemplo `ash 1 g`) conta zero carbono e por isso pode
    esconder um desvio -- a coluna `trocas_sem_formula` avisa."""
    sol = m.optimize()
    if sol.status != "optimal":
        return {"estado": sol.status}
    entra = saiu = 0.0
    sem_formula = []
    for r in m.exchanges:
        f = sol.fluxes.get(r.id, 0.0)
        if abs(f) < 1e-9:
            continue
        met = next(iter(r.metabolites))
        if not met.formula:
            sem_formula.append(r.id)
            continue
        nC = _carbonos(m, r.id)
        if nC == 0:
            continue
        if f < 0:
            entra += -f * nC
        else:
            saiu += f * nC
    mu = sol.fluxes.get(bio_rid, 0.0)
    C_bio = mu * C_por_gDW
    hexose = sum(-sol.fluxes[r.id] for r in m.exchanges
                 if sol.fluxes.get(r.id, 0.0) < 0 and _carbonos(m, r.id) == 6)
    etoh = sol.fluxes.get(etanol_rid, 0.0) if etanol_rid else 0.0
    return {
        "mu": round(mu, 6),
        "cresce": bool(mu > 1e-6),
        "C_entra_mmol": round(entra, 2),
        "C_sai_mmol": round(saiu, 2),
        "C_biomassa_mmol": round(C_bio, 2),
        "fecho_pct": round((saiu + C_bio) / entra * 100, 1) if entra > 1e-9 else None,
        "hexose_consumida": round(hexose, 2),
        "etanol": round(etoh, 2),
        "etanol_por_hexose_mol": round(etoh / hexose, 3) if hexose > 1e-9 else None,
        "pct_do_maximo_teorico": round(etoh / (2 * hexose) * 100, 1) if hexose > 1e-9 else None,
        "trocas_sem_formula": ",".join(sem_formula) or "-",
    }


# ===========================================================================
# T15 -- rendimento em etanol em funcao do mu, e orcamento de carbono
# ===========================================================================
def t15_rendimento(m, bio_rid, etanol_rid, aplicar_meio,
                   mus=(0.0, 0.05, 0.10, 0.20, 0.30, None), q_o2=0.0):
    """Rendimento em etanol com o mu IMPOSTO, e nao no optimo de crescimento.

    Porque e preciso: em anaerobiose a unica fonte de ATP e a glicolise, por
    isso qualquer modelo estequiometrico fermenta tudo o que a estequiometria
    permitir. Ler o rendimento ao mu maximo nao compara capacidade fermentativa,
    compara mu maximos -- e o mu maximo e fixado pela composicao da biomassa e
    pelo GAM, que nestes ficheiros vem do template e nao da especie.

    `mus`: valores de mu a impor. `None` e o mu maximo, para referencia.
    Numa fermentacao vinica o mu anda entre 0.05 e 0.20 h-1 no arranque e cai a
    zero na fase estacionaria; e ai que o rendimento se deve ler.

    Em todas as linhas se fixa o mu e se maximiza o etanol, que e a pergunta
    certa: 'a este ritmo de crescimento, quanto etanol pode este modelo fazer?'.
    A linha `max` fixa o mu no maximo e faz o mesmo, por isso o valor pode ser
    maior do que o etanol que o T13/T14 mostram: la o objectivo e o crescimento
    e o etanol e um entre varios optimos.

    `pct_do_maximo_teorico` acima de 100 nao e erro: o tecto de 2 mol/mol conta
    so a hexose, e o meio de trabalho traz carbono nas vitaminas e nos
    esteroids."""
    linhas = []
    for alvo in mus:
        with m:
            aplicar_meio(m, q_o2=q_o2)
            mumax = _num(m.slim_optimize())
            if alvo is not None and alvo > mumax + 1e-9:
                continue
            mu_fixo = mumax if alvo is None else alvo
            m.reactions.get_by_id(bio_rid).bounds = (mu_fixo, mu_fixo)
            m.objective = etanol_rid
            sol = m.optimize()
            if sol.status != "optimal":
                continue
            f = sol.fluxes
            hexose = sum(-f[r.id] for r in m.exchanges
                         if f.get(r.id, 0.0) < -1e-9 and _carbonos(m, r.id) == 6)
            glic = sum(f[r.id] for r in m.exchanges
                       if (next(iter(r.metabolites)).name or "").lower() == "glycerol")
            etoh = f.get(etanol_rid, 0.0)
            linhas.append({
                "mu_imposto": "max" if alvo is None else alvo,
                "mu": round(f[bio_rid], 4),
                "hexose": round(hexose, 2),
                "etanol": round(etoh, 2),
                "etanol_por_hexose": round(etoh / hexose, 3) if hexose > 1e-9 else None,
                "pct_do_maximo_teorico": round(etoh / (2 * hexose) * 100, 1) if hexose > 1e-9 else None,
                "glicerol": round(glic, 2),
            })
    return pd.DataFrame(linhas)


def t15_orcamento(m, bio_rid, aplicar_meio, q_o2=0.0, limiar=1e-4):
    """Para onde foi cada mmol de carbono, no optimo de crescimento.

    A leitura: a diferenca de rendimento entre dois modelos e igual a diferenca
    no dreno redox (glicerol, poliois, lactato) mais a diferenca no carbono que
    foi para biomassa. `mmolC_por_gDW` sai do proprio balanco e e a coluna que
    diz se duas biomassas sao sequer comparaveis -- um template com menos
    carbono por grama deixa mais carbono para etanol ao mesmo mu, sem que isso
    diga nada sobre a especie."""
    from cobra.flux_analysis import pfba
    with m:
        aplicar_meio(m, q_o2=q_o2)
        sol = pfba(m)
        f, mu = sol.fluxes, sol.fluxes[bio_rid]
        entra = sum(-f[r.id] * _carbonos(m, r.id) for r in m.exchanges
                    if f.get(r.id, 0.0) < -1e-9)
        destinos = {}
        for r in m.exchanges:
            v = f.get(r.id, 0.0)
            if v > limiar and _carbonos(m, r.id) > 0:
                nome = (next(iter(r.metabolites)).name or r.id)[:18]
                destinos[nome] = destinos.get(nome, 0.0) + v * _carbonos(m, r.id)
    saiu = sum(destinos.values())
    destinos["BIOMASSA"] = entra - saiu
    linhas = [{"destino": k, "mmol_C": round(v, 2), "pct": round(v / entra * 100, 1)}
              for k, v in sorted(destinos.items(), key=lambda x: -x[1])]
    resumo = {"mu": round(mu, 4), "C_entra_mmol": round(entra, 2),
              "mmolC_por_gDW": round((entra - saiu) / mu, 1) if mu > 1e-9 else None}
    return resumo, pd.DataFrame(linhas)
