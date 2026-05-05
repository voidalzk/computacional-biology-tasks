library(igraph)

# ============================================================================
# BLOCO INICIAL - Configuracao geral e pastas de trabalho
# ============================================================================

obter_diretorio_script <- function() {
  argumentos <- commandArgs(trailingOnly = FALSE)
  marcador <- "--file="
  linha_arquivo <- argumentos[startsWith(argumentos, marcador)]

  if (length(linha_arquivo) > 0) {
    caminho_script <- sub(marcador, "", linha_arquivo[1])
    return(dirname(normalizePath(caminho_script, winslash = "/", mustWork = FALSE)))
  }

  return(getwd())
}

diretorio_base <- obter_diretorio_script()
arquivo_dataset <- file.path(diretorio_base, "bio-yeast.mtx")
diretorio_figuras <- file.path(diretorio_base, "figuras")

if (!dir.exists(diretorio_figuras)) {
  dir.create(diretorio_figuras, recursive = TRUE)
}

unlink(file.path(diretorio_figuras, "fig_tmp_*.png"), force = TRUE)

ppi <- 200

if (!interactive()) {
  graphics.off()
}

# ============================================================================
# BLOCO AUXILIAR - Funcoes para salvar em imagens PNG
# ============================================================================

arquivo_png_valido <- function(caminho) {
  if (!file.exists(caminho)) {
    return(FALSE)
  }

info <- file.info(caminho)
if (is.na(info$size) || info$size < 8) {
    return(FALSE)
  }

con <- file(caminho, open = "rb")
on.exit(close(con), add = TRUE)
assinatura <- readBin(con, what = "raw", n = 8)

assinatura_ok <- as.raw(c(0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A))
identical(assinatura, assinatura_ok)
}

salvar_png_seguro <- function(caminho_final, plot_fn, largura = 3200, altura = 2200, res = 200) {
  caminho_tmp <- tempfile(pattern = "fig_tmp_", tmpdir = tempdir(), fileext = ".png")

  on.exit({
    if (file.exists(caminho_tmp)) {
      unlink(caminho_tmp, force = TRUE)
    }
  }, add = TRUE)

  png(
    filename = caminho_tmp,
    width = largura,
    height = altura,
    res = res,
    bg = "white"
  )

  tryCatch(
    {
      plot_fn()
    },
    error = function(e) {
      stop(paste("Erro ao desenhar a figura", basename(caminho_final), ":", conditionMessage(e)))
    },
    finally = {
      if (dev.cur() > 1) {
        dev.off()
      }
    }
  )

  if (!arquivo_png_valido(caminho_tmp)) {
    stop(paste("Falha ao gerar PNG valido:", basename(caminho_final)))
  }

  if (file.exists(caminho_final)) {
    suppressWarnings(unlink(caminho_final))
  }

  copiado <- file.copy(caminho_tmp, caminho_final, overwrite = TRUE)
  if (!copiado) {
    caminho_alternativo <- file.path(
      dirname(caminho_final),
      paste0(
        tools::file_path_sans_ext(basename(caminho_final)),
        "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".png"
      )
    )

    copiado_alt <- file.copy(caminho_tmp, caminho_alternativo, overwrite = TRUE)
    if (!copiado_alt) {
      stop(paste("Nao foi possivel salvar a figura:", basename(caminho_final)))
    }

    warning(
      paste(
        "Arquivo alvo bloqueado. Figura salva em:",
        caminho_alternativo
      )
    )
    return(caminho_alternativo)
  }

  caminho_final
}

# ============================================================================
# ETAPA 1. IMPORTACAO E LIMPEZA
# - Importa dataset .mtx
# - Converte para objeto de grafo
# - Identifica se e direcionada
# - Remove loops e arestas multiplas
# ============================================================================

# Leitura do formato Matrix Market (.mtx): remove comentarios (%)
# e descarta a primeira linha util (dimensoes da matriz).
linhas <- readLines(arquivo_dataset, warn = FALSE)
linhas_sem_comentario <- linhas[!startsWith(linhas, "%")]

if (length(linhas_sem_comentario) < 2) {
  stop("Arquivo .mtx invalido: nao ha arestas para carregar.")
}

linhas_arestas <- linhas_sem_comentario[-1]
dados_brutos <- read.table(
  text = linhas_arestas,
  header = FALSE,
  stringsAsFactors = FALSE
)

if (ncol(dados_brutos) < 2) {
  stop("Arquivo .mtx invalido: menos de 2 colunas nas arestas.")
}

# Transformacao para igraph.
network <- as.matrix(dados_brutos[, 1:2])
grafo <- graph_from_edgelist(network, directed = FALSE)

if (is.null(V(grafo)$name)) {
  V(grafo)$name <- as.character(seq_len(vcount(grafo)))
}

grafo <- simplify(grafo, remove.multiple = TRUE, remove.loops = TRUE)

# Identificacao do tipo da rede (direcionada ou nao-direcionada) ainda na Etapa 1.
cat("Direcionado:", is_directed(grafo), "\n")

# ============================================================================
# ETAPA 2. CARACTERIZACAO TOPOLOGICA
# - Ordem, tamanho, densidade, diametro e transitividade global
# ============================================================================

cat("----- REDE COMPLETA -----\n")
cat("Ordem (Vertices):", vcount(grafo), "\n")
cat("Tamanho (Arestas):", ecount(grafo), "\n")
cat("Densidade:", edge_density(grafo, loops = FALSE), "\n")
cat("Diametro:", diameter(grafo, directed = FALSE, unconnected = TRUE), "\n")
cat("Transitividade (Clustering):", transitivity(grafo, type = "globalundirected"), "\n")

# Analise de Grau Basica
dg <- degree(grafo)
cat("Grau medio:", mean(dg), "\n")

# Removendo nós isolados (grau zero)
grafo_sem_isolados <- delete_vertices(grafo, v = degree(grafo) == 0)
cat("\n----- APÓS REMOVER NÓS ISOLADOS -----\n")
cat("Ordem (Vertices):", vcount(grafo_sem_isolados), "\n")
cat("Densidade:", edge_density(grafo_sem_isolados, loops = FALSE), "\n")

# ============================================================================
# ETAPA 3. IDENTIFICACAO DE ATORES-CHAVE (HUBS) E DISTRIBUIÇÃO
# - Centralidade de grau e de intermediacao
# - Ranking Top 5 em cada medida
# ============================================================================

graus <- degree(grafo, mode = "all")
intermediacoes <- betweenness(grafo, directed = FALSE, normalized = FALSE)

top_5_grau_idx <- order(graus, decreasing = TRUE)[1:5]
top_5_inter_idx <- order(intermediacoes, decreasing = TRUE)[1:5]

top_5_grau <- data.frame(
  no = V(grafo)$name[top_5_grau_idx],
  grau = graus[top_5_grau_idx],
  row.names = NULL
)

top_5_inter <- data.frame(
  no = V(grafo)$name[top_5_inter_idx],
  intermediacao = intermediacoes[top_5_inter_idx],
  row.names = NULL
)

cat("\nTop 5 nos por Grau:\n")
print(top_5_grau)

cat("\nTop 5 nos por Intermediacao:\n")
print(top_5_inter)

# ============================================================================
# ETAPA 4. DETECCAO DE COMUNIDADES E VISUALIZACAO
# - Louvain para comunidades
# - Layout Fruchterman-Reingold
# - Tamanho do no proporcional ao grau
# - Cor do no proporcional a comunidade
# ============================================================================

set.seed(42)
comunidades <- cluster_louvain(grafo)
layout_fr <- layout_with_fr(grafo, niter = 1500)

tamanhos_nos <- 2 + 7 * (graus - min(graus)) / (max(graus) - min(graus) + .Machine$double.eps)
membros <- membership(comunidades)
n_comunidades <- length(unique(membros))
paleta <- hcl.colors(n_comunidades, palette = "Dark 3")
cores_nos <- paleta[membros]

arquivo_rede_completa <- file.path(diretorio_figuras, "rede_completa_comunidades_200ppi.png")

arquivo_rede_completa_salvo <- salvar_png_seguro(
  caminho_final = arquivo_rede_completa,
  res = ppi,
  plot_fn = function() {
    par(mar = c(1, 1, 3, 1))
    plot(
      grafo,
      layout = layout_fr,
      vertex.size = tamanhos_nos,
      vertex.color = adjustcolor(cores_nos, alpha.f = 0.9),
      vertex.label = NA,
      vertex.frame.color = NA,
      edge.color = adjustcolor("gray40", alpha.f = 0.12),
      edge.width = 0.4,
      main = paste0("bio-yeast (rede completa) | Louvain: ", n_comunidades, " comunidades")
    )
  }
)

componentes <- components(grafo)
indice_componente_gigante <- which.max(componentes$csize)
vertices_gc <- V(grafo)[componentes$membership == indice_componente_gigante]
grafo_gc <- induced_subgraph(grafo, vids = vertices_gc)

graus_gc <- degree(grafo_gc)
comunidades_gc <- cluster_louvain(grafo_gc)
membros_gc <- membership(comunidades_gc)
n_comunidades_gc <- length(unique(membros_gc))
paleta_gc <- hcl.colors(n_comunidades_gc, palette = "Dark 3")
cores_gc <- paleta_gc[membros_gc]
tamanhos_gc <- 2 + 8 * (graus_gc - min(graus_gc)) / (max(graus_gc) - min(graus_gc) + .Machine$double.eps)

set.seed(42)
layout_gc <- layout_with_fr(grafo_gc, niter = 1800)
hubs_gc_idx <- order(graus_gc, decreasing = TRUE)[1:5]
rotulos_gc <- rep(NA, vcount(grafo_gc))
rotulos_gc[hubs_gc_idx] <- V(grafo_gc)$name[hubs_gc_idx]

arquivo_componente_gigante <- file.path(diretorio_figuras, "componente_gigante_comunidades_200ppi.png")

arquivo_componente_gigante_salvo <- salvar_png_seguro(
  caminho_final = arquivo_componente_gigante,
  res = ppi,
  plot_fn = function() {
    par(mar = c(1, 1, 3, 1))
    plot(
      grafo_gc,
      layout = layout_gc,
      vertex.size = tamanhos_gc,
      vertex.color = adjustcolor(cores_gc, alpha.f = 0.92),
      vertex.label = rotulos_gc,
      vertex.label.cex = 0.7,
      vertex.label.color = "black",
      vertex.frame.color = NA,
      edge.color = adjustcolor("gray35", alpha.f = 0.15),
      edge.width = 0.45,
      main = paste0("bio-yeast (componente gigante) | Louvain: ", n_comunidades_gc, " comunidades")
    )
  }
)

top_10_grau_idx <- order(graus, decreasing = TRUE)[1:10]
top_10_inter_idx <- order(intermediacoes, decreasing = TRUE)[1:10]

arquivo_hubs <- file.path(diretorio_figuras, "top_hubs_grau_intermediacao_200ppi.png")

arquivo_hubs_salvo <- salvar_png_seguro(
  caminho_final = arquivo_hubs,
  largura = 3600,
  altura = 2200,
  res = ppi,
  plot_fn = function() {
    par(mfrow = c(1, 3), mar = c(4, 8, 3, 1), oma = c(0, 0, 1, 0))

    barplot(
      rev(graus[top_10_grau_idx]),
      names.arg = rev(V(grafo)$name[top_10_grau_idx]),
      horiz = TRUE,
      las = 1,
      col = "#2A9D8F",
      border = NA,
      main = "Top 10 por Grau",
      xlab = "Grau"
    )

    barplot(
      rev(intermediacoes[top_10_inter_idx]),
      names.arg = rev(V(grafo)$name[top_10_inter_idx]),
      horiz = TRUE,
      las = 1,
      col = "#E76F51",
      border = NA,
      main = "Top 10 por Intermediacao",
      xlab = "Betweenness"
    )
    
    # Adicionando visualizacao da distribuicao de graus
    hist(dg, breaks = 30, col = "#E9C46A", border = "white",
         main = "Distribuicao de Graus", xlab = "Grau", ylab = "Frequencia")

    mtext("Analise Estrutural e Hubs em bio-yeast", outer = TRUE, cex = 1.1, font = 2)
  }
)

cat("\nImagens salvas em 200 ppi:\n")
cat("-", arquivo_rede_completa_salvo, "\n")
cat("-", arquivo_componente_gigante_salvo, "\n")
cat("-", arquivo_hubs_salvo, "\n")