#' Este script acompanha o "Exercício 1" e realiza o carregamento e o 
#' pré-processamento dos dados indicados na tarefa EaD1, sendo utilizado
#' na etapa seguinte (EaD1 - Script 2).

################################################################################
################################################################################
### SUGESTÃO DE PRÉ-PROCESSAMENTO
################################################################################
################################################################################

################################################################################
### Load data
################################################################################
data_df <- read.table(file = "./WHO-COVID-19-global-table-data.csv", 
    sep=",", quote="\"", header = TRUE, stringsAsFactors = FALSE,
  encoding = "ASCII//TRANSLIT")
colnames(data_df)
# [1] "Name"                                                        
# [2] "WHO.Region"                                                  
# [3] "Cases...cumulative.total"                                    
# [4] "Cases...cumulative.total.per.100000.population"              
# [5] "Cases...newly.reported.in.last.7.days"                       
# [6] "Cases...newly.reported.in.last.7.days.per.100000.population" 
# [7] "Cases...newly.reported.in.last.24.hours"                     
# [8] "Deaths...cumulative.total"                                   
# [9] "Deaths...cumulative.total.per.100000.population"             
# [10] "Deaths...newly.reported.in.last.7.days"                      
# [11] "Deaths...newly.reported.in.last.7.days.per.100000.population"
# [12] "Deaths...newly.reported.in.last.24.hours"

#-------------------------------------------------------------------------------
# Ajuste dos nomes das variáveis
colnames(data_df) <- gsub("...", ".", colnames(data_df), fixed = T)
class(data_df)
"data.frame"

head(data_df)[1:3,1:3]
#  Name      WHO.Region Cases.cumulative.total
# 1       Belarus          Europe                   994084
# 2         China Western Pacific                 99381761
# 3 French Guiana                                    98041

# Verificação das classes das variáveis
class(data_df$Name)
class(data_df$WHO.Region)

#-------------------------------------------------------------------------------
# Preparação de rótulos simplificados para visualização

#--- Rótulos das variáveis
colnames(data_df)
var_labels <- c(
  "Name",
  "Region",
  "Cases_total",
  "Cases_per100k",
  "Cases_7d",
  "Cases_7d_per100k",
  "Cases_24h",
  "Deaths_total",
  "Deaths_per100k",
  "Deaths_7d",
  "Deaths_7d_per100k",
  "Deaths_24h"
)
names(var_labels) <- colnames(data_df)
colnames(data_df) <- var_labels

#--- Rótulos dos registros (linhas)
reg_labels <- substr(data_df$Name, 1, 15)
reg_labels <- make.unique(reg_labels)
names(reg_labels) <- data_df$Name
rownames(data_df) <- reg_labels

################################################################################
### Identificação de variáveis numéricas
################################################################################

#-------------------------------------------------------------------------------
# Verificação das classes com laço 'for'
for(i in 1:ncol(data_df)){
  tp <- class(data_df[,i])
  print(tp)
}
# [1] "character"
# [1] "character"
# [1] "numeric"
# [1] "numeric"
# [1] "numeric"
# [1] "numeric"
# [1] "numeric"
# [1] "numeric"
# [1] "numeric"
# [1] "numeric"
# [1] "logical"
# [1] "numeric"

#-------------------------------------------------------------------------------
# Verificação das classes com 'lapply'
classes_df <- lapply(data_df, class)
classes_df <- unlist(classes_df)

#-------------------------------------------------------------------------------
# Sumário das variáveis numéricas
summary(data_df[, classes_df=="numeric"])

################################################################################
### Salvamento dos dados processados
################################################################################

#-------------------------------------------------------------------------------
# Salvar para análises posteriores
save(data_df, classes_df, var_labels, reg_labels,
  file = "WHO_COVID_19_global_20260323.RData")

