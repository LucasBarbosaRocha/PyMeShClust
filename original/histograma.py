# ###################################################################################
# Classe Histograma
# Esse classe é uma adaptação da classe da biblioteca khmer do python
# Objetivo: Converter uma sequências em um k-mer histogram
# Motivo: Precisei modificar alguns pontos o histograma para ficar o mais próximo do
#         histograma utilizado no artigo, ainda não está tão similar, precisa de alguns 
#         ajustes, mas it is work.
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2010-2015, Michigan State University.
# Copyright (C) 2015-2016, The Regents of the University of California.
# Contact: khmer-project@idyll.org
# pylint: disable=invalid-name,missing-docstring
# ###################################################################################

# ###################################################################################
# importações
# ###################################################################################
import khmer
from khmer import khmer_args
import threading
import numpy as np
import os
from tempfile import NamedTemporaryFile

class Histograma:

    nome = ""
    length = 0
    nkmers = 0
    tablesize = 0
    n_unique_kmers = 0

    def __init__(self, nome, n_unique_kmers, nkmers, tablesize, length):
        self.nome = nome
        self.nkmers = nkmers
        self.tablesize = tablesize
        self.n_unique_kmers = n_unique_kmers
        self.histo = np.zeros((tablesize,4), dtype=np.float32)
        self.length = length

    def imprimeTabela(self):
        print("Sequência " + self.nome)
        print('abundance', 'count', 'cumulative','cumulative_fraction')
        print(self.histo)
        print("Qtd. Máxima k-mers " + str(self.nkmers))
        print("k-mers unicos " + str(self.n_unique_kmers))

# Converte uma sequencia para histograma
def sequenceToHistograma(nome, sequence):

    ksize = 3
    nkmers = 4**ksize
    tablesize = nkmers + 10


    # Initialize countgraph
    cg = khmer.Countgraph(ksize, tablesize, 1)
    # print('Created a countgraph with', cg.hashsizes(), 'buckets')

    # start loading
    # auxNome = "sequenciaAuxliar.fa"
    # aux = open(auxNome, 'w')
    # aux.write(nome)
    # aux.write(sequence+"\n")
    # aux.close()

    # fp = TemporaryFile('w+t')
    # fp = ff.TemporaryFile(mode='w+t', suffix=".fasta")

    with NamedTemporaryFile(prefix="lucas", suffix=".fasta", delete=False, mode="w+t") as fp:
        fp.write(nome)
        fp.write(sequence)
        fp.seek(0)

    # print(fp.name)

    # print(fp.read())
    # print(fp.read())
    fp.close()
    rparser = khmer.ReadParser(fp.name)
    # fp.close()
    # os.remove(fp.name)
    # os.unlink(fp.name)

    # rparser2 = rparser
    # aux.closes
    # os.remove(auxNome)
    # rparser = khmer.ReadParser(sequence)
    threads = []
    for _ in range(1):
        thread = \
            threading.Thread(
                target=cg.consume_seqfile_with_reads_parser,
                args=(rparser, )
            )
        threads.append(thread)
        thread.start()
    
    for thread in threads:
        thread.join()

    # print('unique', cg.n_unique_kmers())
    h = Histograma(nome, cg.n_unique_kmers(), nkmers, tablesize, len(sequence))

    abundance_lists = []
    
    tracking = khmer_args.create_matching_nodegraph(cg)
    def __do_abundance_dist__(read_parser):
        abundances = cg.abundance_distribution_with_reads_parser(
            read_parser, tracking)
        abundance_lists.append(abundances)

    # with NamedTemporaryFile(prefix="lucas", suffix=".fasta", delete=False, mode="w+t") as fp:
    #     fp.write(nome)
    #     fp.write(sequence)
    #     fp.seek(0)

    # print(fp.name)

    rparser2 = khmer.ReadParser(fp.name)
    # fp.close()
    # os.remove(fp.name)
    # # os.unlink(fp.name)

    threads = []

    for _ in range(1):
        thread = \
            threading.Thread(
                target=__do_abundance_dist__,
                args=(rparser2, )
            )
        threads.append(thread)
        thread.start()

    for thread in threads:
        thread.join()

    assert len(abundance_lists) == 1, len(abundance_lists)
    abundance = {}
    for abundance_list in abundance_lists:
        for i, count in enumerate(abundance_list):
            abundance[i] = abundance.get(i, 0) + count

    total = sum(abundance.values())

    if 0 == total:
        print("ERROR: abundance distribution is uniformly zero; "
                  "nothing to report.")
        print("\tPlease verify that the input files are valid.")
        return 0

    sofar = 0
    line = 0
    for _, i in sorted(abundance.items()):
        if i == 0 and line < h.tablesize:
            continue

        sofar += i
        frac = sofar / float(total)

        #hist_fp_csv.writerow([_, i, sofar, round(frac, 3)])
        # print(line, tablesize, [_, i, sofar, round(frac, 3)])
        h.histo[line][0] = _
        h.histo[line][1] = i
        h.histo[line][2] = sofar
        h.histo[line][3] = round(frac, 3)
        line = line + 1
        if sofar == total:
            break
    
    return h