
#' Este script acompanha o "Exercício 2", carregando no RStudio sugestões de  
#' dados indicados na tarefa EaD2, juntamente com sugestões de análise.
#' 
#' Parte da rede "ca-Erdos992" para realizar algumas análises estruturais e 
#' visualizações comparativas. Use como um ponto de partida para suas análises
#' e lembre de complementar com as demais informações descritas na tarefa EaD2.

library(dplyr)
library(ggplot2)
library(stringr)
library(plyr)
library(RColorBrewer)
library(igraph)
library(RGraphSpace)

################################################################################
### Carregamento da tabela de dados "ca-Erdos992"
### https://networkrepository.com/ca-Erdos992.php
################################################################################
Erdos992  <- read.table("ca-Erdos992/ca-Erdos992.mtx", 
    header = F, comment.char = "%")

head(Erdos992)
#    V1 V2
# 1 183  1
# 2 299  1
# 3 309  1
# 4 384  1
# 5 411  1
# 6 432  1

################################################################################
### Conversão dos dados em um objeto igraph
################################################################################
network <- as.matrix(Erdos992 [,1:2])
network <- igraph::graph_from_edgelist(network, directed = FALSE)
class(network)
# [1] "igraph"

#-------------------------------------------------------------------------------
# Verifica a ordem e o tamanho (número total de vértices e arestas)
vcount(network)
# [1] 6100
ecount(network)
# [1] 7515

#-------------------------------------------------------------------------------
# Análise do grau dos nós (nodos)
dg <- degree(network)
#--- Grau médio
mean(dg)
# [1] 2.463934
#--- Os 5 graus mais frequentes
table(sort(dg))[1:5]
#    0    1    2    3    4 
# 1006 3600  633  229  121
#--- Histograma de distribuição de graus
hist(dg)

#-------------------------------------------------------------------------------
# Verifica a densidade (proporção entre conexões existentes e possíveis)
edge_density(network)
# [1] 0.0004039899

#-------------------------------------------------------------------------------
# Verifica o diâmetro (maior distância geodésica entre dois nós)
diameter(network)
# [1] 14

#-------------------------------------------------------------------------------
# Verifica a transitividade (coeficiente de agrupamento/clustering coefficient)
# Define a probabilidade de que os vértices adjacentes de um nó estejam conectados
transitivity(network, type = "global")
# [1] 0.04194019

################################################################################
### Visualização do grafo com RGraphSpace
################################################################################

#-------------------------------------------------------------------------------
# Ajusta parâmetros estéticos da rede
V(network)$nodeSize <- 2
E(network)$edgeLineWidth <- 0.3
E(network)$edgeLineColor <- "lightblue"

#-------------------------------------------------------------------------------
# Gera a plotagem do grafo
gs1 <- GraphSpace(network, layout = layout_with_drl(network))
plotGraphSpace(gs1)

#-------------------------------------------------------------------------------
# Remove nós de grau zero (isolados) e refaz o grafo
network2 <- delete_vertices(network, v = degree(network)==0 )

gs2 <- GraphSpace(network2, layout = layout_with_drl(network2))
plotGraphSpace(gs2)

#-------------------------------------------------------------------------------
# Extrai o subgrafo contendo o maior componente conexo
network3 <- largest_component(network)

gs3 <- GraphSpace(network3, layout = layout_with_drl(network3))
plotGraphSpace(gs3)

################################################################################
### Reavaliação das estatísticas (Pós-filtragem)
################################################################################

#-------------------------------------------------------------------------------
# Verifica nova ordem e tamanho
vcount(network3)
# [1] 4991 #era 6100
ecount(network3)
# [1] 7428 #era 7515

#-------------------------------------------------------------------------------
# Reavalia o grau dos nós
dg <- degree(network3)
#--- Novo grau médio
mean(dg)
# [1] 2.976558 #era 2.463934
#--- Novos 5 graus mais frequentes
table(sort(dg))[1:5]
#    1    2    3    4    5 
# 3514  629  225  117   67 
#--- Novo histograma
hist(dg)

#-------------------------------------------------------------------------------
# Verifica a nova densidade
edge_density(network3)
# [1] 0.0005965046 #era 0.0004039899

#-------------------------------------------------------------------------------
# Verifica o diâmetro (permanece o mesmo se o caminho mais longo estiver no componente principal)
diameter(network3)
# [1] 14 #era 14

#-------------------------------------------------------------------------------
# Verifica a nova transitividade
transitivity(network3, type = "global")
# [1] 0.04204424 #era 0.04194019

#-------------------------------------------------------------------------------
#' Ao comparar os resultados, observamos que a remoção de 1.006 nodos 
#' isolados (grau zero) e a filtragem pelo maior componente conexo
#' tornaram a rede mais densa e conectada. Embora o número total de 
#' vértices e arestas tenha diminuído, o grau médio subiu de 2,46 para 
#' 2,97, indicando que os nodos restantes são, em média, mais conectados. 
#' Curiosamente, o diâmetro permaneceu em 14, o que revela que a maior 
#' distância geodésica da rede original já estava contida em seu núcleo 
#' principal, enquanto a transitividade teve um leve aumento, sugerindo 
#' uma estrutura de agrupamento um pouco mais coesa após a limpeza dos dados.

