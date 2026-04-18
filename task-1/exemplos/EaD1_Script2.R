#' Este script acompanha o "Exercício 1" e realiza o carregamento dos 
#' dados previamente pré-processados na etapa anterior (EaD1 - Script 1).

################################################################################
################################################################################
### Análise Exploratória de Dados de COVID-19 (WHO)
################################################################################
################################################################################

################################################################################
### Preparação dos dados para visualização em mapas de calor
################################################################################
library(ComplexHeatmap)
library(RColorBrewer)

#-------------------------------------------------------------------------------
# Carrega dados pré-processados
load("WHO_COVID_19_global_20260323.RData")

#-------------------------------------------------------------------------------
# Seleciona subconjunto de interesse (exemplo: Europa + variáveis numéricas)
data_df2 <- data_df[data_df$Region=="Europe" , classes_df=="numeric"]

#-------------------------------------------------------------------------------
# Verifica a quantidade de valores ausentes (NAs) por variável
# Observe que algumas variáveis apresentam muitos NAs
apply(is.na(data_df2), 2, sum)
# Cases.cumulative.total 
# 0 
# Cases.cumulative.total.per.100000.population 
# 9 
# Cases.newly.reported.in.last.7.days 
# 180 

#-------------------------------------------------------------------------------
# Remove variáveis com grande quantidade de NAs
data_df2 <- data_df2[, apply(is.na(data_df2), 2, sum) < 20]

#-------------------------------------------------------------------------------
# Mantém apenas observações completas (sem NAs)
data_df2 <- data_df2[complete.cases(data_df2), ]

dim(data_df2)
# [1] 37  6

#-------------------------------------------------------------------------------
# Padroniza os dados (z-score) e gera mapa de calor com agrupamento hierárquico
X <- scale(data_df2)

Heatmap(X, 
  row_names_gp = gpar(fontsize = 7), 
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  clustering_distance_columns = "euclidean",
  clustering_distance_rows = "euclidean")

################################################################################
### Preparação dos dados para visualização em grafo de associação
### (conforme exemplo da tarefa EaD1)
################################################################################
library(RedeR)
library(TreeAndLeaf)
library(RColorBrewer)
library(igraph)

#-------------------------------------------------------------------------------
# Clusteriza os registros (agrupamento de países)
# Etapas: padronização das variáveis e cálculo da distância entre registros
data_scaled <- scale(data_df2)

hc <- hclust(dist(data_scaled, method = "euclidean"), "ward.D2")
plot(hc, cex=0.3)

#-------------------------------------------------------------------------------
# Converte o objeto 'hclust' para objeto 'treeAndLeaf'
tal <- treeAndLeaf(hc)

#-------------------------------------------------------------------------------
# Mapeia variáveis de interesse no objeto 'tal'
tal <- att.mapv(g = tal, dat = data_df2, refcol = 0)

pal <- rev(brewer.pal(7, "RdYlGn"))

tal <- att.setv(g = tal, from = "Deaths_total", 
  to = "nodeColor", cols = pal, nquant = 7)

tal <- att.setv(g = tal, from = "Cases_total",
  to = "nodeSize", xlim = c(10, 50, 5), nquant = 7)

#-------------------------------------------------------------------------------
# Ajusta o tamanho da fonte dos nós folha
V(tal)$nodeFontSize[V(tal)$isLeaf] <- 300

#-------------------------------------------------------------------------------
# Inicializa a interface de visualização
startRedeR()
resetRedeR()

#-------------------------------------------------------------------------------
# Adiciona o grafo à interface
addGraphToRedeR(g = tal, zoom = 75)

# Seleciona nós internos (não folhas) como âncoras
selectNodes(V(tal)$name[!V(tal)$isLeaf], anchor=TRUE)

# Adiciona legendas
addLegendToRedeR(tal)
addLegendToRedeR(tal, type = "nodesize")

# Ajusta o layout (relaxamento do grafo)
relaxRedeR(p1 = 10, p3 = 1, p4 = 100, p5 = 1)

