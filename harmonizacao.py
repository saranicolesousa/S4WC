"""
harmonizacao.py — reconciliacao do compartimento extracelular entre o yeast-GEM
e os tres modelos CarveFungi, antes de qualquer analise de comunidade.

Porque e necessario
-------------------
O SMETANA constroi o meio partilhado agrupando as trocas dos varios organismos
pelo identificador do metabolito extracelular. Os quatro modelos deste projeto
usam dois espacos de nomes incompativeis:

    yeast-GEM    metabolito s_0565     troca r_1714     ("D-glucose")
    CarveFungi   metabolito glc__D_e   troca UF02xxx_E  ("D-glucose")

Sem reconciliacao, a glucose da S. cerevisiae e a glucose da T. delbrueckii sao
dois pools distintos. O consorcio fica particionado: nenhuma troca entre a
Saccharomyces e as nao-Saccharomyces e detetavel, e todos os indices (MIP, MRO,
SC, MUS, MPS) sao calculados sobre uma comunidade que nao partilha meio nenhum.
Os resultados nao ficam errados por pouco: ficam estruturalmente vazios,
exatamente na classe de interacao que este projeto quer medir.

Estrategia
----------
O espaco de nomes de destino e o dos CarveFungi (estilo BiGG, `<id>_e`), por
duas razoes: tres dos quatro modelos ja o usam, e e o que o `reframed` espera
com `flavor='bigg'`.

A correspondencia usa tres chaves independentes, por ordem de prioridade:

    1. bigg.metabolite   anotacao do yeast-GEM + sufixo `_e`
    2. ChEBI             normalizado (o yeast-GEM escreve `CHEBI:`, os
                         CarveFungi escrevem `ChEBI:`)
    3. KEGG compound

Nao se usa o nome como chave primaria. Os nomes divergem em casos que importam
("D-glucose" vs "glucose", "L-argininium(1+)" vs "L-arginine") e coincidem em
casos que nao deviam. O nome entra apenas no `CURADAS`, revisto a mao.

O ponto que distingue esta implementacao de uma tabela de traducao cega: as
tres chaves sao calculadas todas, e os desacordos entre elas sao devolvidos em
`conflitos` em vez de serem resolvidos em silencio pela prioridade. Um conflito
e um erro de anotacao num dos dois modelos e tem de ser lido, nao arrumado.
"""

from __future__ import annotations
import re
from collections import defaultdict


# --------------------------------------------------------------------------
# normalizacao de chaves
# --------------------------------------------------------------------------

def _norm_chebi(v):
    """`CHEBI:16004`, `ChEBI:16004`, `16004` -> `16004`."""
    if v is None:
        return None
    return str(v).strip().upper().replace("CHEBI:", "")


def _norm_kegg(v):
    if v is None:
        return None
    return str(v).strip().upper()


def _primeiro(v):
    """As anotacoes COBRApy sao str ou list[str]. Devolve o primeiro valor."""
    if v is None:
        return None
    if isinstance(v, (list, tuple, set)):
        v = sorted(v)
        return v[0] if v else None
    return v


def _todos(v):
    if v is None:
        return []
    if isinstance(v, (list, tuple, set)):
        return [str(x) for x in v]
    return [str(v)]


# --------------------------------------------------------------------------
# correspondencias revistas a mao
# --------------------------------------------------------------------------
# So entram aqui compostos em que as tres chaves automaticas falham e a
# identidade quimica e inequivoca. Cada linha tem de ter uma razao escrita.

CURADAS = {
    # yeast-GEM -> CarveFungi    razao da correspondencia manual
    #
    # Os tres casos abaixo tem identidade quimica inequivoca (mesma formula,
    # mesmo composto) e falham as tres chaves automaticas porque o metabolito
    # CarveFungi so traz anotacao PubChem, que o yeast-GEM nao usa.
    "s_1295": "9Zhxda_e",   # palmitoleato = (9Z)-hexadecenoato, C16H29O2,
                            # PubChem 445638. Auxotrofia identificada no obj2.
    "s_2826": "9Zocda_e",   # oleato = (9Z)-octadec-9-enoato, C18H33O2,
                            # PubChem 5460221. E o suplemento do Tween 80.
    "s_4242": "3mppal_e",   # metional = 3-(metiltio)propanal, C4H8OS.
                            # Produto de Ehrlich a partir da metionina.
}

# Correspondencias plausiveis que ficaram deliberadamente por fazer, com a
# razao. Entram no relatorio do notebook para nao passarem por esquecimento.
POR_VERIFICAR = {
    "s_4249": ("4hpheetoh_e",
               "tirosol = 4-(2-hidroxietil)fenol; o metabolito CarveFungi "
               "chama-se 4-(1-hidroxietil)fenol. Mesma formula, posicao do "
               "hidroxilo diferente. Um dos dois nomes esta errado no ficheiro "
               "e nao ha anotacao que decida. Por confirmar na reacao que o "
               "produz antes de juntar os pools."),
    "s_4045": ("dtmp_e",
               "timidina 3'-monofosfato; o yeast-GEM tem tambem o 5' (s_4047) "
               "e os CarveFungi so tem um dTMP. Juntar exige decidir qual, e "
               "nenhum dos dois e relevante para troca em mosto."),
}

# Metabolitos extracelulares que NAO devem ser reconciliados mesmo que uma
# chave sugira correspondencia. Sao pseudo-metabolitos de ficheiro, nao
# especies quimicas, e junta-los cria um pool partilhado sem sentido fisico.
NUNCA = {
    "ash_1g_e",         # cinza mineral agregada dos CarveFungi (sem formula)
}


# --------------------------------------------------------------------------
# construcao do mapa
# --------------------------------------------------------------------------

def _indexar_alvo(modelos_alvo):
    """Indexa os metabolitos extracelulares dos modelos de destino pelas tres
    chaves. Devolve (por_bigg, por_chebi, por_kegg, universo)."""
    por_bigg, por_chebi, por_kegg = {}, defaultdict(set), defaultdict(set)
    universo = {}
    for m in modelos_alvo:
        for met in m.metabolites:
            if met.compartment != "e" or met.id in NUNCA:
                continue
            universo[met.id] = met
            por_bigg[met.id] = met.id            # o proprio id ja e BiGG+_e
            for c in _todos(met.annotation.get("chebi")):
                por_chebi[_norm_chebi(c)].add(met.id)
            for c in _todos(met.annotation.get("kegg.compound")):
                por_kegg[_norm_kegg(c)].add(met.id)
    return por_bigg, por_chebi, por_kegg, universo


def construir_mapa_extracelular(modelo_fonte, modelos_alvo, verbose=True):
    """Constroi a correspondencia s_#### -> <bigg>_e.

    Devolve um dicionario com:
      mapa       {id_fonte: id_alvo}                 correspondencias aceites
      chave      {id_fonte: 'bigg'|'chebi'|'kegg'|'curada'}
      conflitos  [(id_fonte, {chave: id_alvo, ...})] chaves em desacordo
      sem_par    [id_fonte, ...]                     sem correspondencia
      orfaos     [id_alvo, ...]                      so existem no destino
    """
    por_bigg, por_chebi, por_kegg, universo = _indexar_alvo(modelos_alvo)

    mapa, chave, conflitos, sem_par = {}, {}, [], []

    for met in modelo_fonte.metabolites:
        if met.compartment != "e":
            continue

        cand = {}

        if met.id in CURADAS:
            cand["curada"] = CURADAS[met.id]

        b = _primeiro(met.annotation.get("bigg.metabolite"))
        if b and f"{b}_e" in por_bigg:
            cand["bigg"] = f"{b}_e"

        for c in _todos(met.annotation.get("chebi")):
            hits = por_chebi.get(_norm_chebi(c))
            if hits and len(hits) == 1:
                cand["chebi"] = next(iter(hits))
                break

        for c in _todos(met.annotation.get("kegg.compound")):
            hits = por_kegg.get(_norm_kegg(c))
            if hits and len(hits) == 1:
                cand["kegg"] = next(iter(hits))
                break

        if not cand:
            sem_par.append(met.id)
            continue

        distintos = set(cand.values())
        if len(distintos) > 1:
            conflitos.append((met.id, dict(cand)))

        for k in ("curada", "bigg", "chebi", "kegg"):
            if k in cand:
                mapa[met.id] = cand[k]
                chave[met.id] = k
                break

    usados = set(mapa.values())
    orfaos = sorted(set(universo) - usados)

    r = {"mapa": mapa, "chave": chave, "conflitos": conflitos,
         "sem_par": sorted(sem_par), "orfaos": orfaos,
         "universo_alvo": sorted(universo)}

    if verbose:
        n_fonte = sum(1 for x in modelo_fonte.metabolites if x.compartment == "e")
        print(f"fonte, metabolitos extracelulares : {n_fonte}")
        print(f"destino, universo                 : {len(universo)}")
        print(f"correspondidos                    : {len(mapa)}")
        for k in ("curada", "bigg", "chebi", "kegg"):
            n = sum(1 for v in chave.values() if v == k)
            if n:
                print(f"   por {k:8s}                   : {n}")
        print(f"sem par (so na fonte)             : {len(sem_par)}")
        print(f"orfaos  (so no destino)           : {len(orfaos)}")
        print(f"conflitos entre chaves            : {len(conflitos)}")

    return r


# --------------------------------------------------------------------------
# aplicacao
# --------------------------------------------------------------------------

_SUF = re.compile(r"_e$")


def aplicar_mapa(modelo, mapa, prefixo_troca="EX_", verbose=False):
    """Renomeia os metabolitos extracelulares do modelo segundo `mapa` e
    reescreve os identificadores das trocas para `EX_<metabolito>`.

    Trabalha sobre o modelo recebido. Passe uma copia se precisar do original.
    Devolve o numero de (metabolitos, trocas) renomeados.
    """
    n_met = 0
    for antigo, novo in mapa.items():
        if antigo not in modelo.metabolites:
            continue
        met = modelo.metabolites.get_by_id(antigo)
        if novo in modelo.metabolites and novo != antigo:
            # colisao: o destino ja existe neste modelo. Nao fundir em silencio.
            raise ValueError(
                f"colisao ao renomear {antigo} -> {novo}: o destino ja existe")
        met.id = novo
        n_met += 1

    # Os extracelulares que ficaram sem par mantem o id de origem, e no
    # yeast-GEM esse id nao traz o compartimento (`s_1374`). O SMETANA recupera
    # o nome do composto com `original_metabolite[2:-2]`, ou seja tira `M_` da
    # frente e DOIS caracteres do fim, assumindo o sufixo `_e`. Sem esse sufixo
    # a etiqueta sai truncada e, pior, varios metabolitos colapsam na mesma:
    # s_4199, s_4200 e s_4201 sairiam todos como `s_42`. Os relatorios de meio
    # minimo ficariam ilegiveis e com compostos fundidos. Acrescenta-se o
    # sufixo para que a convencao do SMETANA se aplique sem perda.
    n_suf = 0
    for met in modelo.metabolites:
        if met.compartment == "e" and not met.id.endswith("_e"):
            alvo = f"{met.id}_e"
            if alvo in modelo.metabolites:
                raise ValueError(f"colisao ao sufixar {met.id} -> {alvo}")
            met.id = alvo
            n_suf += 1
    modelo.repair()

    n_tr = 0
    for r in list(modelo.exchanges):
        mets = list(r.metabolites)
        if len(mets) != 1:
            continue
        novo_id = f"{prefixo_troca}{mets[0].id}"
        if novo_id == r.id:
            continue
        if novo_id in modelo.reactions:
            raise ValueError(f"colisao ao renomear troca {r.id} -> {novo_id}")
        r.id = novo_id
        n_tr += 1
    modelo.repair()

    if verbose:
        print(f"   {n_met} metabolitos renomeados, {n_suf} sufixados com _e, "
              f"{n_tr} trocas renomeadas")
    return n_met, n_tr


def normalizar_trocas(modelo, prefixo_troca="EX_", verbose=False):
    """So reescreve os ids das trocas para `EX_<metabolito>`, sem mexer nos
    metabolitos. Usa-se nos modelos que ja estao no espaco de nomes de destino
    (os CarveFungi, cujas trocas se chamam `UF#####_E`)."""
    return aplicar_mapa(modelo, {}, prefixo_troca=prefixo_troca, verbose=verbose)


def diagnostico_pools(modelos, nomes=None):
    """Quantos metabolitos extracelulares sao partilhados por cada par de
    modelos. E a medida directa de quanto meio comum o SMETANA vai ver."""
    import itertools
    nomes = nomes or [m.id for m in modelos]
    ext = {n: {x.id for x in m.metabolites if x.compartment == "e"}
           for n, m in zip(nomes, modelos)}
    linhas = []
    for a, b in itertools.combinations(nomes, 2):
        inter = ext[a] & ext[b]
        linhas.append({"A": a, "B": b, "ext_A": len(ext[a]), "ext_B": len(ext[b]),
                       "partilhados": len(inter),
                       "jaccard": len(inter) / len(ext[a] | ext[b])})
    return linhas
