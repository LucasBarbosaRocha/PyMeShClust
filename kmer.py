# ###################################################################################
# Módulo kmer
# Objetivo: ler arquivos FASTA e converter cada sequência em um vetor de frequência
#           de k-mers (4^k posições: AAA, AAC, ..., TTT para k=3).
# Memória:  todas as sequências ficam em uma única matriz numpy (N, 4^k) float32,
#           pré-alocada. Para k=3 são 256 bytes por sequência.
# ###################################################################################

import numpy as np

# Tabela byte -> código da base (A=0, C=1, G=2, T=3, qualquer outro = -1)
_CODIGO = np.full(256, -1, dtype=np.int64)
for _i, _b in enumerate(b"ACGT"):
    _CODIGO[_b] = _i
    _CODIGO[_b + 32] = _i  # minúsculas


def ler_fasta(caminho):
    """Gera (nome, sequencia) sem carregar o arquivo inteiro.

    Aceita sequências quebradas em várias linhas. O nome vem sem o '>'.
    """
    nome, partes = None, []
    with open(caminho) as arquivo:
        for linha in arquivo:
            linha = linha.strip()
            if not linha:
                continue
            if linha.startswith(">"):
                if nome is not None:
                    yield nome, "".join(partes)
                nome, partes = linha[1:], []
            else:
                partes.append(linha)
    if nome is not None:
        yield nome, "".join(partes)


def contar_sequencias(caminho):
    with open(caminho) as arquivo:
        return sum(1 for linha in arquivo if linha.startswith(">"))


def vetor_kmer(sequencia, k=3):
    """Conta os k-mers da sequência. Janelas com base diferente de ACGT são ignoradas."""
    codigos = _CODIGO[np.frombuffer(sequencia.encode("ascii", "replace"), dtype=np.uint8)]
    n = len(codigos) - k + 1
    if n <= 0:
        return np.zeros(4**k, dtype=np.int64)

    indice = np.zeros(n, dtype=np.int64)
    valido = np.ones(n, dtype=bool)
    for j in range(k):
        janela = codigos[j:j + n]
        valido &= janela >= 0
        indice = indice * 4 + janela
    return np.bincount(indice[valido], minlength=4**k)


def carregar(caminho, k=3):
    """Lê o FASTA em duas passadas (conta, depois preenche a matriz pré-alocada).

    Retorna (nomes, comprimentos, X), onde cada linha de X é a frequência
    normalizada dos k-mers (soma 1). Normalizar remove o viés do comprimento.
    Sequências sem nenhum k-mer válido ficam com vetor zero.
    """
    n = contar_sequencias(caminho)
    X = np.zeros((n, 4**k), dtype=np.float32)
    comprimentos = np.zeros(n, dtype=np.int64)
    nomes = []
    for i, (nome, sequencia) in enumerate(ler_fasta(caminho)):
        contagem = vetor_kmer(sequencia, k)
        total = contagem.sum()
        if total:
            X[i] = contagem / total
        comprimentos[i] = len(sequencia)
        nomes.append(nome)
    return nomes, comprimentos, X


def agrupar(labels):
    """Devolve (grupos, ids): grupos[i] são os índices das sequências do cluster ids[i]."""
    ordem = np.argsort(labels, kind="stable")
    cortes = np.flatnonzero(np.diff(labels[ordem])) + 1
    ids = labels[ordem][np.r_[0, cortes]] if len(ordem) else np.array([], dtype=labels.dtype)
    return np.split(ordem, cortes), ids


def escrever_clstr(caminho, nomes, comprimentos, labels, representantes):
    """Escreve no formato .clstr do CD-HIT. O representante do cluster recebe '*'.

    representantes[c] é o índice da sequência representante do cluster c.
    """
    grupos, ids = agrupar(labels)
    with open(caminho, "w") as saida:
        for grupo, c in zip(grupos, ids):
            saida.write(f">Cluster {c}\n")
            for j, i in enumerate(grupo):
                marca = " *" if i == representantes[c] else ""
                saida.write(f"{j}\t{comprimentos[i]}nt, >{nomes[i]}{marca}\n")
