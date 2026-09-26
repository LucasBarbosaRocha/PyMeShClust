# ###################################################################################
# Experimentos do PyMeShClust
# Roda o clustering variando o tamanho do k-mer (3 a 6) e o limiar de similaridade,
# avalia contra os rótulos dos cabeçalhos e escreve RESULTADOS.md e resultados.csv.
#
# Uso: python experimentos.py
# ###################################################################################

import csv
import os
import time
import tracemalloc

from avaliar import ROTULOS, avaliar, ler_clstr
from kmer import carregar
from meshclust import meshclust

KMERS = [3, 4, 5, 6]
LIMIARES = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
DATASETS = [
    # arquivo, rótulos avaliados e se agrupa só as sequências com rótulo
    # (remove os reads curtos sem rótulo)
    ("hbv.fasta", ["genotipo", "regiao"], False),
    ("sequencias.fasta", ["especie", "genero"], False),
    ("sequencias.fasta", ["especie", "genero"], True),
]
ORIGINAL = [
    ("sequencias.fasta", "original/output95.clstr", "original (2020)", 0.95),
    ("sequencias.fasta", "original/output80.clstr", "original (2020)", 0.80),
]


def rodar(arquivo, kmer, limiar, filtro=None):
    tracemalloc.start()
    inicio = time.perf_counter()
    nomes, _, X = carregar(arquivo, kmer)
    if filtro:
        manter = [filtro(nome) is not None for nome in nomes]
        nomes = [nome for nome, m in zip(nomes, manter) if m]
        X = X[manter]
    labels = meshclust(X, limiar)[0]
    tempo = time.perf_counter() - inicio
    pico = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    return nomes, labels, tempo, pico


def formatar(m):
    return f"{m['ari']:.3f} / {m['nmi']:.3f} / {m['pureza']:.3f}"


def linha_csv(titulo, metodo, kmer, limiar, rotulo, m, tempo="", pico=""):
    return [titulo, metodo, kmer, limiar, rotulo, m["clusters"], m["rotuladas"], m["classes"],
            f"{m['ari']:.4f}", f"{m['nmi']:.4f}", f"{m['pureza']:.4f}", tempo, pico]


def main():
    linhas_csv = []
    md = ["# Resultados do PyMeShClust", "",
          "Gerado por `python experimentos.py`. Cada célula de rótulo mostra **ARI / NMI / pureza** "
          "(1 = perfeito; ARI ~0 = aleatório). Só as sequências com rótulo no cabeçalho entram na "
          "avaliação, mas todas são agrupadas. Memória = pico de alocação do numpy/Python "
          "(tracemalloc) para ler o FASTA e rodar o clustering (com a fase de junção).", ""]

    for arquivo, rotulos, so_rotuladas in DATASETS:
        filtro = ROTULOS[rotulos[0]] if so_rotuladas else None
        titulo = f"{arquivo} (só sequências com rótulo)" if so_rotuladas else arquivo
        md += [f"## {titulo}", "",
               "| k-mer | limiar | clusters | " + " | ".join(rotulos) + " | tempo (s) | memória (MB) |",
               "|---|---|---|" + "---|" * len(rotulos) + "---|---|"]
        for kmer in KMERS:
            for limiar in LIMIARES:
                nomes, labels, tempo, pico = rodar(arquivo, kmer, limiar, filtro)
                metricas = {r: avaliar(nomes, labels, r) for r in rotulos}
                clusters = metricas[rotulos[0]]["clusters"]
                md.append(f"| {kmer} | {limiar} | {clusters} | "
                          + " | ".join(formatar(metricas[r]) for r in rotulos)
                          + f" | {tempo:.2f} | {pico / 2**20:.1f} |")
                for r in rotulos:
                    linhas_csv.append(linha_csv(titulo, "meshclust", kmer, limiar, r, metricas[r],
                                                f"{tempo:.3f}", f"{pico / 2**20:.2f}"))
                print(f"{titulo} k-mer={kmer} limiar={limiar} ({tempo:.2f}s)")

        for arq, clstr, descricao, limiar in ORIGINAL:
            if arq == arquivo and not so_rotuladas and os.path.exists(clstr):
                nomes, clusters = ler_clstr(clstr)
                metricas = {r: avaliar(nomes, clusters, r) for r in rotulos}
                md.append(f"| {descricao} | {limiar} | {metricas[rotulos[0]]['clusters']} | "
                          + " | ".join(formatar(metricas[r]) for r in rotulos) + " | - | - |")
                for r in rotulos:
                    linhas_csv.append(linha_csv(titulo, "original", "", limiar, r, metricas[r]))
        md.append("")

    with open("RESULTADOS.md", "w") as saida:
        saida.write("\n".join(md))
    with open("resultados.csv", "w", newline="") as saida:
        escritor = csv.writer(saida)
        escritor.writerow(["arquivo", "metodo", "kmer", "limiar", "rotulo", "clusters", "rotuladas",
                           "classes", "ari", "nmi", "pureza", "tempo_s", "memoria_mb"])
        escritor.writerows(linhas_csv)
    print("==> RESULTADOS.md e resultados.csv criados!")


if __name__ == "__main__":
    main()
