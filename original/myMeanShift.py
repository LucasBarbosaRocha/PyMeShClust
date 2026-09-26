# ###################################################################################
# Classe MeanShfit
# Essa classe é uma implementação do MeanShift utilizando sklearn
# Objetivo: Dados vários histrogramas, achar o melhor centroid dos histogramas e então 
#           trocar o centroid atual para o novo centroid
# ###################################################################################

# ###################################################################################
# importações
# ###################################################################################
from cluster import *
from sklearn.cluster import MeanShift, estimate_bandwidth
import numpy as np
import cv2 as cv
import math
from scipy.stats import pearsonr

# Calcula a similaridade entre dois histogramas
def return_intersection(hist_1, hist_2):
    minima = np.minimum(hist_1, hist_2)
    intersection = np.true_divide(np.sum(minima), np.sum(hist_2))
    return intersection

def return_intersection2(hist_1, hist_2):
    minima = np.minimum(hist_1, hist_2) * 2
    intersection = np.true_divide(np.sum(minima), np.sum(hist_2 + hist_1))
    return intersection

def mySigmoid(x):
    return 1/(1 + np.exp(-x))

def calculandoSimilaride(length_A, length_B, sequencia_A, sequencia_B, k, similaridade):
    myC = return_intersection(sequencia_A[0],sequencia_B[0])
    if (myC >= similaridade):
        return True
    return False

# Calcula a similaridade entre uma sequência e um cluster e devolve a sequência mais similar
def devolveNovoCentroidComIntersection(centroid_antigo_cluster, sequencias, similaridade):
    centroid_aux = sequencias[0]
    old_taxa = 0
    for i in sequencias:
        if (i != centroid_aux):
            # myP = myPearson(centroid_antigo_cluster.centroid.histo, i.histo)
            myP = return_intersection(centroid_antigo_cluster.centroid.histo, i.histo)
            if (myP > old_taxa):
                old_taxa = myP
                centroid_aux = i
    return centroid_aux

# Executando o meanshift para criar um centroid sintetico e então calcular
# a sequência mais próxima do centroid sintetitco 
# e trocar o centroid atual
def executaMeanShift2(cluster_atual, similaridade, kmer):
    # print("#### MeanShift Rodando...")
    x = []
    # print ("Cluster atual contém: ", len(cluster_atual.cluster))
    for k in cluster_atual.cluster:
        if (k != cluster_atual.centroid):
            aux = np.array(k.histo.T[1:2][0])
            for i in range(len(aux)):
                if ([i, aux[i]] not in x):
                    x.append([i,aux[i]])
    
    aux = cluster_atual.centroid.histo.T[1:2][0]
    p = []
    for i in range(len(aux)):
        p.append([i,aux[i]])

    ms = MeanShift(bandwidth=2, bin_seeding=True)
    ms.fit(x)
    predit = np.array([ms.predict(p)])

    centroid_novo = cluster_atual.centroid
    x = cluster_atual.centroid.histo.T[1:2][0]
    p = predit[0]
    aux = []
    for i in range(len(x)):
        if (x[i] == 0):
            aux.append(0)
        else:
            aux.append(p[i])

    maior = return_intersection(np.array(aux)[0], cluster_atual.centroid.histo.T[1:2][0])
    for k in cluster_atual.cluster:
        x = k.histo.T[1:2][0]

        aux = []
        for i in range(len(x)):
            if (x[i] == 0):
                aux.append(0)
            else:
                aux.append(p[i])
        x = k.histo.T[1:2]

        r = return_intersection(np.array(aux)[0], x[0])
        if (r > maior):
            centroid_novo = k # Atualizando o centroid
            maior =r
    return centroid_novo
