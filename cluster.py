# #################################################################
# Classe Cluster
# Essa classe é uma representação do cluster utilizado no trabalho.
# Objetivo: Construir um cluster de histrogramas, onde cada objeto
#           objeto cluster possui um centroid e vários histogramas
#           que pertencem ao atual cluster
# #################################################################

# ###################################################################################
# importações
# ###################################################################################
from histograma import *

class Cluster:

    def __init__(self, centroid, sequencia):
        self.centroid = centroid
        self.cluster = []
        self.cluster.append(sequencia)

    def insereCluster(self, sequencia):
        self.cluster.append(sequencia)   
    
    def trocaCentroid(self, centroid):
        self.centroid = centroid
        if not centroid in self.cluster:
            self.cluster.append(centroid)
            