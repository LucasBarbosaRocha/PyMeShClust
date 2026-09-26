# ###################################################################################
# Nome: Lucas Barbosa Rocha
# Disciplina: Inteligência Artificial
# Trabalho: Implementar um clustering para sequências de DNA utilizando MeanShift.
# Contato: lucas.lb.rocha@gmail.com
# Git: Lucasbarbosarocha
#
# Objetivo: implementar o trabalho "MeShClust: an intelligent tool for clustering DNA
#           sequences". Cada sequência vira um vetor de frequência de k-mers; os
#           clusters são formados por mean shift guloso e depois juntados.
#           A saída segue o formato .clstr do CD-HIT.
#
# Uso: python main.py sequencias.fasta -s 0.95
# ###################################################################################

import argparse
import time

from kmer import carregar, escrever_clstr
from meshclust import meshclust


def main():
    parser = argparse.ArgumentParser(description="Clustering de sequências de DNA estilo MeShClust.")
    parser.add_argument("entrada", nargs="?", default="sequencias.fasta", help="arquivo FASTA")
    parser.add_argument("-s", "--similaridade", type=float, default=0.95,
                        help="limiar de similaridade de k-mers entre 0 e 1")
    parser.add_argument("--kmer", type=int, default=3, help="tamanho do k-mer (memória: 4^k floats por sequência)")
    parser.add_argument("--sem-juncao", action="store_true", help="pula a fase de junção de clusters")
    parser.add_argument("-o", "--saida", help="arquivo .clstr (padrão: output<similaridade>.clstr)")
    args = parser.parse_args()
    if not 0 < args.similaridade <= 1:
        parser.error("similaridade deve estar entre 0 e 1")
    saida = args.saida or f"output{round(100 * args.similaridade)}.clstr"

    inicio = time.perf_counter()
    print("### Convertendo sequências para vetores de k-mers.")
    nomes, comprimentos, X = carregar(args.entrada, args.kmer)
    print(f"==> {len(nomes)} sequências convertidas ({X.nbytes / 2**20:.2f} MB).")

    print("### Gerando clusters.")
    labels, centros = meshclust(X, args.similaridade, juntar=not args.sem_juncao)
    print(f"==> {len(centros)} clusters gerados.")

    print("### Escrevendo no arquivo de saída.")
    escrever_clstr(saida, nomes, comprimentos, labels, centros)
    print(f"==> {saida} criado! ({time.perf_counter() - inicio:.2f}s)")


if __name__ == "__main__":
    main()
