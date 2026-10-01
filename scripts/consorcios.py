"""
consorcios.py — construcao das comunidades e traducao do MS300 para o ambiente partilhado do SMETANA.
-------------------------
O `meio_mosto.py` aplica o MS300 a um modelo isolado, escrevendo o bound de cada especie na troca correspondente. 
Numa comunidade do SMETANA ha dois niveis de restricao:

    pool   `R_EX_M_<met>_pool`      o que o MEIO fornece ao consorcio inteiro
    troca  `R_EX_<met>_<org>`       o que CADA organismo pode tirar do pool

Dar so o pool nao é suficiente: sem limite por organismo, um dos modelos pode consumir
a totalidade de um nutriente.
Em termos de oxigénio, o regime não é o mesmo para todos: a S. cerevisiae, a T. delbrueckii e a L. thermotolerans
correm em anaerobiose (q_O2 = 0) e a M. pulcherrima em microaerobiose
(q_O2 = 0.5). O pool abre O2 ao valor da M. pulcherrima e as trocas individuais dos outros tres ficam a zero. 
-> testar com 0.5 para todos posteriormente
"""

from __future__ import annotations
import itertools
import meio_mosto as MM


# ── nomes ────────────────────────────────────────────────────────────────────

CODIGOS = {"Sc": "Scer", "Td": "Tdel", "Mp": "Mpul", "Lt": "Ltho"}
ESPECIES = {
    "Sc": "Saccharomyces cerevisiae",
    "Td": "Torulaspora delbrueckii",
    "Mp": "Metschnikowia pulcherrima",
    "Lt": "Lachancea thermotolerans",
}

Q_O2 = {"Sc": 0.0, "Td": 0.0, "Mp": 0.5, "Lt": 0.0}

# ── mapa especie do MS300 -> id do metabolito harmonizado ───────────────────

def mapa_especie_metabolito(modelos_cobra):
    """{especie do MS300: id do metabolito extracelular} apos harmonizacao.

    `modelos_cobra` e {codigo: modelo cobrapy JA harmonizado}. Percorre-se
    todos porque nem todas as especies existem em todos os modelos, e verifica-se
    que os que a tem concordam no id.
    """
    mapa, divergencias = {}, {}
    for cod, m in modelos_cobra.items():
        mm = MM.construir_mapa(m)
        for especie, rid in mm.items():
            if rid is None:
                continue
            met = next(iter(m.reactions.get_by_id(rid).metabolites)).id
            if especie in mapa and mapa[especie] != met:
                divergencias.setdefault(especie, {}).update(
                    {cod: met, "anterior": mapa[especie]})
            mapa.setdefault(especie, met)
    for cod, m in modelos_cobra.items():
        for r in m.exchanges:
            met = next(iter(r.metabolites))
            if MM._norm(met.name) in ("oxygen", "dioxygen", "o2"):
                mapa.setdefault("oxygen", met.id)
                break

    return mapa, divergencias


# ── ambiente ─────────────────────────────────────────────────────────────────
def _pool(met_id):
    return f"R_EX_M_{met_id}_pool"


def _troca(met_id, org_id):
    return f"R_EX_{met_id}_{org_id}"


def ambiente_mosto(comunidade, mapa_met, codigos=None, mu=MM.MU_REF,
                   x_final=MM.X_REF, eficiencia=None, suplementos=None,
                   q_suplemento=None, ash=None, verbose=False):
    """Constroi o dicionario de bounds do pool para o MS300.
    """
    from reframed import Environment

    suplementos = MM.SUPLEMENTOS if suplementos is None else suplementos
    q_suplemento = MM.Q_SUPLEMENTO if q_suplemento is None else q_suplemento
    bounds = MM.bounds_do_meio(mu=mu, x_final=x_final, eficiencia=eficiencia)
    ash_b = (MM.calcular_bound(MM.ash_mmol_L(), mu, x_final)
             if ash is None else -abs(ash))

    merged = comunidade.merged
    amb, fora, aplicados = {}, [], {}

    def por(met_id, lb):
        """Um composto entra no ambiente se o bound for negativo.
        """
        p = _pool(met_id)
        if p not in merged.reactions or abs(lb) < 1e-12:
            return False
        amb[p] = (lb, 1000.0)
        return True

    # camada 2 - basais sem limite
    for k in MM.BASAIS:
        met = mapa_met.get(k)
        if met and por(met, -1000.0):
            aplicados[k] = (met, -1000.0)
        else:
            fora.append((k, "basal sem metabolito no merged"))

    # camada 2b - suplementos declarados, limite medido
    for k in suplementos:
        met = mapa_met.get(k)
        if met and por(met, -abs(q_suplemento)):
            aplicados[k] = (met, -abs(q_suplemento))
        else:
            fora.append((k, "suplemento sem metabolito no merged"))

    # camada 3 - cinza mineral dos CarveFungi
    met = mapa_met.get("ash")
    if met and por(met, ash_b):
        aplicados["ash"] = (met, ash_b)

    # camada 4 - composicao do mosto
    for especie, b in bounds.items():
        if especie in MM.BASAIS or especie in suplementos:
            continue
        met = mapa_met.get(especie)
        if met is None:
            fora.append((especie, "sem metabolito apos harmonizacao"))
            continue
        if not por(met, b):
            fora.append((especie, f"{met} sem pool no merged"))
            continue
        aplicados[especie] = (met, b)

    # camada 5 - oxigenio
    met = mapa_met.get("oxygen") or mapa_met.get("dioxygen")

    q_max = max(Q_O2[c] for c in (codigos or Q_O2))
    if met and por(met, -abs(q_max)):
        aplicados["oxygen"] = (met, -abs(q_max))

    env = Environment()
    env.update(amb)

    rel = {"aplicados": aplicados, "fora": fora,
           "n_pool": len(amb), "ash": ash_b, "q_o2_pool": -abs(q_max)}
    if verbose:
        print(f"  pools abertos: {len(amb)}")
        print(f"  fora: {len(fora)} -> {[f[0] for f in fora]}")
    return env, rel


def preparar_para_comunidade(modelo, cod, q_o2=None, abertura=1000.0):
    """Poe um modelo ja harmonizado no estado em que entra numa comunidade."""
    q_o2 = Q_O2[cod] if q_o2 is None else q_o2
    o2, n = None, 0
    for r in modelo.exchanges:
        nome = MM._norm(next(iter(r.metabolites)).name)
        if nome in ("oxygen", "dioxygen", "o2"):
            o2 = r.id
            continue
        r.lower_bound = -abs(abertura)
        n += 1
    if o2:
        modelo.reactions.get_by_id(o2).lower_bound = -abs(q_o2)
    return {"trocas_abertas": n, "o2": o2, "q_o2": -abs(q_o2)}


# ── desenhos de consorcio ────────────────────────────────────────────────────

def desenhos(codigos=("Sc", "Td", "Mp", "Lt"), tamanhos=(2, 4)):
    """consorcios a testar.
    Pares e o quarteto. -> testar com trios posteriormente
    """
    out = {}
    for t in tamanhos:
        for combo in itertools.combinations(codigos, t):
            out["+".join(combo)] = list(combo)
    return out
