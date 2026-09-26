# ###################################################################################
# Módulo avaliar
# Objetivo: medir a qualidade dos clusters comparando com um rótulo verdadeiro
#           extraído do cabeçalho de cada sequência.
#
# Rótulos:
#   genotipo  genótipo do HBV ("HBV genotype C ..."  -> C; F2 -> F)
#   regiao    região do genoma do HBV (genoma completo, gene S, gene P, P/preS/S)
#   especie   gênero + espécie ("Klebsiella pneumoniae strain ..." -> Klebsiella pneumoniae)
#   genero    só o gênero
# Sequências sem rótulo (ex.: reads "ERR599000.3-...") ficam fora da avaliação.
#
# Métricas (1 = perfeito):
#   ARI     Adjusted Rand Index: concordância de pares, corrigida pelo acaso (~0 = aleatório)
#   NMI     informação mútua normalizada entre clusters e classes
#   pureza  fração das sequências que pertencem à classe majoritária do seu cluster
#           (sobe artificialmente com muitos clusters pequenos; olhar junto com o ARI)
#
# Uso: python avaliar.py output95.clstr --rotulo especie
# ###################################################################################

import argparse
import re

import numpy as np


def rotulo_genotipo(nome):
    m = re.search(r"genotype ([A-H])", nome)
    return m.group(1) if m else None


def rotulo_regiao(nome):
    if "complete genome" in nome or "complete cds" in nome:
        return "genoma completo"
    if "preS1" in nome:
        return "P/preS/S parcial"
    if " S gene" in nome:
        return "gene S parcial"
    if "Pol gene" in nome or "polymerase (P) gene" in nome:
        return "gene P parcial"
    return None


def rotulo_especie(nome):
    partes = nome.split()
    if len(partes) >= 3 and re.fullmatch(r"[A-Z][a-z]+", partes[1]) and partes[2].islower():
        return f"{partes[1]} {partes[2]}"
    return None


def rotulo_genero(nome):
    especie = rotulo_especie(nome)
    return especie.split()[0] if especie else None


ROTULOS = {
    "genotipo": rotulo_genotipo,
    "regiao": rotulo_regiao,
    "especie": rotulo_especie,
    "genero": rotulo_genero,
}


def ler_clstr(caminho):
    """Lê um .clstr (formato CD-HIT) e devolve (nomes, clusters) na ordem do arquivo."""
    nomes, clusters = [], []
    cluster = -1
    with open(caminho) as arquivo:
        for linha in arquivo:
            if linha.startswith(">Cluster"):
                cluster += 1
            elif ">" in linha:
                nome = linha.split(">", 1)[1].rstrip("\n")
                if nome.endswith(" *"):
                    nome = nome[:-2]
                nomes.append(nome)
                clusters.append(cluster)
    return nomes, np.array(clusters)


def contingencia(classes, clusters):
    _, c = np.unique(classes, return_inverse=True)
    _, k = np.unique(clusters, return_inverse=True)
    tabela = np.zeros((c.max() + 1, k.max() + 1), dtype=np.int64)
    np.add.at(tabela, (c, k), 1)
    return tabela


def _pares(x):
    return x * (x - 1) / 2


def ari(tabela):
    n = tabela.sum()
    soma_ij = _pares(tabela).sum()
    soma_a = _pares(tabela.sum(axis=1)).sum()
    soma_b = _pares(tabela.sum(axis=0)).sum()
    esperado = soma_a * soma_b / _pares(n) if n > 1 else 0.0
    maximo = (soma_a + soma_b) / 2
    if maximo == esperado:
        return 1.0
    return (soma_ij - esperado) / (maximo - esperado)


def _entropia(contagens):
    p = contagens[contagens > 0] / contagens.sum()
    return -(p * np.log(p)).sum()


def nmi(tabela):
    n = tabela.sum()
    a, b = tabela.sum(axis=1), tabela.sum(axis=0)
    i, j = np.nonzero(tabela)
    nij = tabela[i, j]
    info_mutua = (nij / n * np.log(nij * n / (a[i] * b[j]))).sum()
    media = (_entropia(a) + _entropia(b)) / 2
    return 1.0 if media == 0 else info_mutua / media


def pureza(tabela):
    return tabela.max(axis=0).sum() / tabela.sum()


def avaliar(nomes, clusters, rotulo):
    """Métricas dos clusters contra o rótulo; só as sequências com rótulo entram."""
    funcao = ROTULOS[rotulo]
    classes = [funcao(nome) for nome in nomes]
    tem = np.array([c is not None for c in classes])
    if not tem.any():
        raise ValueError(f"nenhuma sequência tem rótulo '{rotulo}'")
    tabela = contingencia(np.array([c for c in classes if c is not None]), np.asarray(clusters)[tem])
    return {
        "rotuladas": int(tem.sum()),
        "classes": tabela.shape[0],
        "clusters": len(np.unique(clusters)),
        "ari": ari(tabela),
        "nmi": nmi(tabela),
        "pureza": pureza(tabela),
    }


def main():
    parser = argparse.ArgumentParser(description="Avalia um .clstr contra rótulos tirados do cabeçalho.")
    parser.add_argument("clstr", help="arquivo .clstr")
    parser.add_argument("-r", "--rotulo", choices=ROTULOS, required=True)
    args = parser.parse_args()

    nomes, clusters = ler_clstr(args.clstr)
    try:
        m = avaliar(nomes, clusters, args.rotulo)
    except ValueError as erro:
        parser.error(str(erro))
    print(f"{m['clusters']} clusters | {m['rotuladas']} de {len(nomes)} sequências com rótulo, "
          f"{m['classes']} classes ({args.rotulo})")
    print(f"ARI {m['ari']:.3f} | NMI {m['nmi']:.3f} | pureza {m['pureza']:.3f}")


if __name__ == "__main__":
    main()
