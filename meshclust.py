# ###################################################################################
# Módulo meshclust
# Objetivo: clustering no estilo MeShClust sobre os vetores de k-mers.
#
# Fase 1 (mean shift guloso):
#   - Começa com um centro. Todas as sequências livres com similaridade >= limiar
#     entram no cluster.
#   - Mean shift: a média dos vetores dos membros é o centro sintético; o novo centro
#     é o membro real mais parecido com essa média. Repete a busca com o novo centro
#     até nenhuma sequência nova entrar.
#   - O próximo cluster começa pela sequência livre mais parecida com o último
#     centro (como no MeShClust).
# Fase 2 (junção):
#   - Clusters cujos centros têm similaridade >= limiar são juntados (do maior para
#     o menor, sem encadear) e o centro é recalculado. Repete até estabilizar.
#
# Memória: tudo trabalha sobre a matriz X (N, 4^k) e um vetor de rótulos; nenhuma
# lista de objetos por sequência. As similaridades são calculadas em blocos para
# não criar temporários do tamanho de X.
# ###################################################################################

import numpy as np

from kmer import agrupar

BLOCO = 65536


def similaridade(A, b):
    """Interseção de histogramas normalizados: soma(min(a, b)). Simétrica, entre 0 e 1."""
    return np.minimum(A, b).sum(axis=-1)


def similaridades(X, indices, b):
    """Similaridade de X[indices] com o vetor b, calculada em blocos."""
    saida = np.empty(len(indices), dtype=np.float32)
    for inicio in range(0, len(indices), BLOCO):
        fim = inicio + BLOCO
        saida[inicio:fim] = similaridade(X[indices[inicio:fim]], b)
    return saida


def centro_do_cluster(X, membros):
    """Membro real mais parecido com a média dos membros (centro sintético)."""
    media = X[membros].mean(axis=0)
    return membros[similaridades(X, membros, media).argmax()]


def fase_mean_shift(X, limiar, inicio=0):
    n = len(X)
    labels = np.full(n, -1, dtype=np.int64)
    centros = []
    livres = np.arange(n)
    semente = inicio

    while livres.size:
        c = len(centros)
        labels[semente] = c
        livres = livres[livres != semente]
        membros = [np.array([semente])]
        soma = X[semente].astype(np.float64)
        quantidade = 1
        centro = semente

        while livres.size:
            sim = similaridades(X, livres, X[centro])
            entra = sim >= limiar
            if not entra.any():
                break
            novos = livres[entra]
            livres = livres[~entra]
            labels[novos] = c
            membros.append(novos)
            soma += X[novos].sum(axis=0, dtype=np.float64)
            quantidade += len(novos)

            todos = np.concatenate(membros)
            membros = [todos]
            media = (soma / quantidade).astype(X.dtype)
            centro = todos[similaridades(X, todos, media).argmax()]

        centros.append(centro)
        if livres.size:
            semente = livres[similaridades(X, livres, X[centro]).argmax()]

    return labels, np.array(centros, dtype=np.int64)


def fase_juncao(X, labels, limiar, max_rodadas=20):
    for _ in range(max_rodadas):
        grupos, _ = agrupar(labels)
        centros = np.array([centro_do_cluster(X, g) for g in grupos])
        destino = np.arange(len(grupos))
        mantidos = []
        for i in sorted(range(len(grupos)), key=lambda i: -len(grupos[i])):
            if mantidos:
                sim = similaridades(X, centros[mantidos], X[centros[i]])
                j = sim.argmax()
                if sim[j] >= limiar:
                    destino[i] = mantidos[j]
                    continue
            mantidos.append(i)

        if len(mantidos) == len(grupos):
            break
        novos = np.empty_like(labels)
        for i, g in enumerate(grupos):
            novos[g] = destino[i]
        labels = novos

    # Renumera 0..C-1 mantendo a ordem de criação dos clusters
    _, labels = np.unique(labels, return_inverse=True)
    grupos, _ = agrupar(labels)
    centros = np.array([centro_do_cluster(X, g) for g in grupos], dtype=np.int64)
    return labels, centros


def meshclust(X, limiar=0.95, juntar=True):
    """Retorna (labels, centros): centros[c] é o índice da sequência centro do cluster c."""
    labels, centros = fase_mean_shift(X, limiar)
    if juntar:
        labels, centros = fase_juncao(X, labels, limiar)
    return labels, centros
