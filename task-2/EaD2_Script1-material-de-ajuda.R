
#' Este script acompanha o "Exercício 2", carregando no RStudio os dados 
#' indicados na tarefa EaD2.

library(dplyr)
library(ggplot2)
library(stringr)
library(plyr)
library(RColorBrewer)
library(igraph)
library(RGraphSpace)

################################################################################
### Carrega a tabela de dados "aves-wildbird-network"
################################################################################
# https://networkrepository.com/
wildbird <- read.table("aves-wildbird-network/aves-wildbird-network.edges", 
    header = F)

################################################################################
### Transforma dados em um objeto igraph
################################################################################
network <- as.matrix(wildbird[,1:2])
network <- igraph::graph_from_edgelist(network, directed = FALSE)
class(network)
# [1] "igraph"
     
gs <- GraphSpace(network, layout = layout_with_kk(network))
plotGraphSpace(gs)


