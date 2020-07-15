# ###################################################################################
# Nome: Lucas Barbosa Rocha
# Disciplina: Inteligência Artificial
# Trabalho: Implementar um clustering para sequências de DNA 
#           utilizando MeanShift.
# Contato: lucas.lb.rocha@gmail.com
# Git: Lucasbarbosarocha
#
# Classe main
# Objetivo: implementar o trabalho: MeShClust: an intelligent tool for clustering DNA
#           sequences. O foco foi implementar (fazer funcionar). Os autores utilizaram
#           c++ e fizeram implementações na mão de muitos pontos. Eu trouxe para o python
#           e tentei utilizar as bibliotecas do python. Alguns pontos não estão 
#           apresentando os mesmos resultados porque ainda estou em processo de
#           interpretação de como os autores estão fazendo no trabalho deles.
# Problema: Cuidado com o tamanho do arquivo de entrada. A utilização da memória não 
#           está otimizada, ou seja, estou criando várias variáveis e o consumo de 
#           memória está alta.
# ###################################################################################

# ###################################################################################
# importações
# ###################################################################################
import os
import requests
import numpy as np
from histograma import *
from abundance_dist_single import *
from cluster import *

from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn import linear_model
from sklearn.linear_model import LogisticRegression
from myMeanShift import *
from sklearn.metrics import accuracy_score, jaccard_score
from sklearn.linear_model import SGDClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import Perceptron
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import TweedieRegressor
from sklearn.linear_model import Ridge
from sklearn.linear_model import ElasticNet
from sklearn.linear_model import GammaRegressor
from sklearn.linear_model import Lasso, LassoCV, LinearRegression
from scipy.stats import poisson
from glimpy import GLM, Poisson
# from sklearn.model_selection import cross_val_scor
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
# from sklearn import datasets
from sklearn.model_selection import GridSearchCV

# ###################################################################################
# main()
# ###################################################################################
# executar()
# Abrindo arquivo
nomeArquivo = "sequencias.fasta"
arquivo = open(nomeArquivo)
similaridade = 0.95
auxNome = "output"+str(int(100*similaridade))+".clstr"
kmer = 3
# Convertendo 'todas' as sequências para histograma
sequencias = []
print ("### Convertendo sequências para k-mer histograma.")
qtd_sequencias = 0
while True and qtd_sequencias < 500:
    nome = arquivo.readline()
    if len(nome) == 0:
        break
    sequencia = arquivo.readline()
    # print (nome)
    histo = sequenceToHistograma(nome, sequencia)
    # histo.imprimeTabela()
    sequencias.append(histo)
    qtd_sequencias = qtd_sequencias + 1
# Fechando arquivo
# arquivo.close()
print ("==> ",qtd_sequencias, " Sequências convertidas!")

# Rodando o clusterizador
i = 0
clusters = []
cluster_atual = Cluster(sequencias[0], sequencias[0])
Centros = []
Clusteres = []
sequencias.remove(sequencias[0])
qtd_c = 0
# Primeira Parte, gerando os clusteres
# Lucas Lembrete: cuidando com a lista não esvaziando
print("### Gerando Clusteres.")
val = 0
g = []
while True: 
    # Procedimento em construção para rodar de 500 em 500 sequências
    if len(sequencias) == 0 and val == 0:
        nome = arquivo.readline()
        if (len(nome) != 0): # Conferindo o resto do arquivo
            print ("### Convertendo sequências para k-mer histograma.")
            qtd_sequencias = 0
            aux = 1
            while True and qtd_sequencias < 500:
                if (aux != 1):
                    nome = arquivo.readline()
                aux = 2
                if len(nome) == 0:
                    break
                sequencia = arquivo.readline()
                # print (nome)
                histo = sequenceToHistograma(nome, sequencia)
                # histo.imprimeTabela()
                sequencias.append(histo)
                qtd_sequencias = qtd_sequencias + 1
            # Fechando arquivo
            # arquivo.close()
            print ("==> ",qtd_sequencias, " Sequência convertidas!")
        else:
            break
    g = []
    i = i + 1
    # Percorrer as sequências e determinar as sequências similares ao centroid atual
    for j in sequencias:
        if (calculandoSimilaride(cluster_atual.centroid.length, j.length, cluster_atual.centroid.histo.T[1:2], j.histo.T[1:2], kmer, similaridade)):
            g.append(j) # sequencias similares com o centroid

    # Para cada sequências similar, colocar no cluster do centroid atual
    for j in g:
        if j not in cluster_atual.cluster: # Jutando e enviando para o cluster atual
            cluster_atual.insereCluster(j)
        sequencias.remove(j) # Removendo sequências de S

    if len(g) != 0:
        # Caso g não esteja vazio, precisamos atualizar o centroid atual
        # rodar meanshift , determinar um centroid sintetico e então escolher um novo centroid
        # utilizando a similaridade de intersecao
        cluster_atual.centroid = executaMeanShift2(cluster_atual, similaridade, kmer)
    else:
        # Caso g fique vazio, um novo cluster deve ser criado
        Centros.append(cluster_atual.centroid)
        Clusteres.append(cluster_atual)
        val = 0
        if len(sequencias) > 0:
            qtd_c = qtd_c + 1
            sequencia_aux = devolveNovoCentroidComIntersection(cluster_atual, sequencias, similaridade)
            cluster_atual = Cluster(sequencia_aux, sequencia_aux)
            sequencias.remove(sequencia_aux)
            val = 1

# conferir se na última iteração o g não esta cheio e saiu porque o S esvaziou
if len(g) != 0:
    Centros.append(cluster_atual.centroid)
    Clusteres.append(cluster_atual)
    if len(sequencias) > 0:
        cluster_atual = Cluster(sequencias[0], sequencias[0])
        sequencias.remove(sequencias[0])   
# print("centros ", Centros)
print("==> Clusteres gerados!")

# imprimindo resultados
print("### Escrevendo no arquivo de saída.")
arquivoSaida = open(auxNome, 'w')
qtd_cluster = 0
for i in Clusteres:
    # print ('>Cluster ', qtd_cluster)
    e = ">Cluster " + str(qtd_cluster)
    arquivoSaida.write(e+"\n")
    aux = 0
    for j in i.cluster:
        if (j.nome == i.centroid.nome):
            # print("%d %dnt, %s *\n" %(aux, j.length, j.nome[0:len(j.nome)-1]), end='') 
            e = str(aux) + " " + str(j.length)+"nt, " + j.nome[0:len(j.nome)-1] + " *"
            arquivoSaida.write(str(e)+"\n")
        else:
            # print("%d %dnt, %s" %(aux, j.length, j.nome), end='') 
            e = str(aux) + " " + str(j.length) +"nt, "+ j.nome[0:len(j.nome)-1]
            arquivoSaida.write(str(e)+"\n")
        aux = aux + 1
    qtd_cluster = qtd_cluster + 1
print("==> Done!")
print ("==>", auxNome, " criado!")
arquivoSaida.close()
