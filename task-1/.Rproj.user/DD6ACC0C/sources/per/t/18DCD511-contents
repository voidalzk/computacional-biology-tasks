# Exercício 1 - Visualizando dados da OMS com abordagem de redes

###############################################################################
# 1) Importação do arquivo CSV e organização em data.frame
###############################################################################

suppressPackageStartupMessages({
	library(igraph)
})

rescale_to <- function(x, to = c(0, 1)) {
	rx <- range(x, na.rm = TRUE)
	if (is.infinite(rx[1]) || is.infinite(rx[2]) || rx[1] == rx[2]) {
		return(rep(mean(to), length(x)))
	}
	(x - rx[1]) / (rx[2] - rx[1]) * (to[2] - to[1]) + to[1]
}

csv_path <- "./WHO-COVID-19-global-table-data.csv"

who_df <- read.csv(
	file = csv_path,
	stringsAsFactors = FALSE,
	check.names = FALSE,
	fileEncoding = "UTF-8-BOM"
)

stopifnot(is.data.frame(who_df))

# Renomeia colunas para facilitar manipulação
col_map <- c(
	"Name" = "country",
	"WHO Region" = "region",
	"Cases - cumulative total" = "cases_total",
	"Cases - cumulative total per 100000 population" = "cases_per100k",
	"Cases - newly reported in last 7 days" = "cases_7d",
	"Cases - newly reported in last 7 days per 100000 population" = "cases_7d_per100k",
	"Cases - newly reported in last 24 hours" = "cases_24h",
	"Deaths - cumulative total" = "deaths_total",
	"Deaths - cumulative total per 100000 population" = "deaths_per100k",
	"Deaths - newly reported in last 7 days" = "deaths_7d",
	"Deaths - newly reported in last 7 days per 100000 population" = "deaths_7d_per100k",
	"Deaths - newly reported in last 24 hours" = "deaths_24h"
)

names(who_df) <- ifelse(names(who_df) %in% names(col_map),
	col_map[names(who_df)],
	names(who_df)
)

num_cols <- c(
	"cases_total", "cases_per100k", "cases_7d", "cases_7d_per100k",
	"cases_24h", "deaths_total", "deaths_per100k", "deaths_7d",
	"deaths_7d_per100k", "deaths_24h"
)

for (cn in num_cols) {
	if (cn %in% names(who_df)) {
		who_df[[cn]] <- suppressWarnings(as.numeric(who_df[[cn]]))
	}
}

who_df$region[is.na(who_df$region) | who_df$region == ""] <- "Unknown"

# Remove linhas agregadas que distorcem comparações entre países
who_df <- subset(
	who_df,
	country != "Global" & !grepl("^International conveyance", country)
)

cat("Dimensão do data.frame:", dim(who_df)[1], "linhas x", dim(who_df)[2], "colunas\n")

###############################################################################
# 2) Estratégias de visualização em rede
###############################################################################

if (!dir.exists("figuras")) dir.create("figuras")

# -----------------------------------------------------------------------------
# Estratégia A: Rede bipartida Região -> País (peso por casos acumulados)
# -----------------------------------------------------------------------------

# Mantém os 8 países com mais casos por região para reduzir poluição visual
split_by_region <- split(who_df, who_df$region)
top_country_df <- do.call(
	rbind,
	lapply(split_by_region, function(d) {
		d <- d[order(d$cases_total, decreasing = TRUE), ]
		head(d, 8)
	})
)

edges_a <- data.frame(
	from = top_country_df$region,
	to = top_country_df$country,
	weight = log10(top_country_df$cases_total + 1),
	cases_total = top_country_df$cases_total,
	stringsAsFactors = FALSE
)

g_a <- graph_from_data_frame(edges_a, directed = FALSE)

V(g_a)$type <- V(g_a)$name %in% unique(edges_a$to)
V(g_a)$color <- ifelse(V(g_a)$type, "#1f77b4", "#d62728")
V(g_a)$size <- ifelse(V(g_a)$type, 7, 16)
V(g_a)$label.cex <- ifelse(V(g_a)$type, 0.45, 0.8)
E(g_a)$width <- rescale_to(E(g_a)$weight, to = c(0.7, 5))
E(g_a)$color <- "#bdbdbd"


png("figuras/rede_bipartida_regiao_pais.png", width = 3200, height = 2200, res = 300)
plot(
	g_a,
	layout = layout_with_fr(g_a),
	main = "Estratégia A: Rede bipartida (Região - País)\nAresta ponderada por log10(casos acumulados + 1)",
	vertex.label.family = "sans",
	vertex.label.color = "black",
	vertex.frame.color = NA
)
legend(
	"topleft",
	legend = c("País", "Região OMS"),
	pt.bg = c("#1f77b4", "#d62728"),
	pch = 21,
	pt.cex = c(1.2, 1.6),
	bty = "n",
	title = "Cores dos nós"
)
dev.off()

# -----------------------------------------------------------------------------
# Estratégia B: Rede de similaridade entre regiões
# -----------------------------------------------------------------------------

region_df <- aggregate(
	who_df[, c("cases_total", "cases_per100k", "deaths_total", "deaths_per100k")],
	by = list(region = who_df$region),
	FUN = function(x) mean(x, na.rm = TRUE)
)

row.names(region_df) <- region_df$region
region_mat <- as.matrix(region_df[, -1])
region_mat <- scale(region_mat)

dist_mat <- as.matrix(dist(region_mat, method = "euclidean"))
sim_mat <- exp(-dist_mat)

diag(sim_mat) <- 0
thr <- quantile(sim_mat[upper.tri(sim_mat)], probs = 0.75, na.rm = TRUE)

edge_idx <- which(sim_mat >= thr, arr.ind = TRUE)
edge_idx <- edge_idx[edge_idx[, 1] < edge_idx[, 2], , drop = FALSE]

if (nrow(edge_idx) == 0) {
	# Fallback: conecta cada região ao vizinho mais parecido
	edge_idx <- do.call(rbind, lapply(seq_len(nrow(sim_mat)), function(i) {
		j <- which.max(sim_mat[i, ])
		c(i, j)
	}))
	edge_idx <- edge_idx[edge_idx[, 1] < edge_idx[, 2], , drop = FALSE]
}

edges_b <- data.frame(
	from = row.names(sim_mat)[edge_idx[, 1]],
	to = row.names(sim_mat)[edge_idx[, 2]],
	similarity = sim_mat[edge_idx],
	stringsAsFactors = FALSE
)

g_b <- graph_from_data_frame(edges_b, directed = FALSE, vertices = region_df$region)

region_total_cases <- aggregate(cases_total ~ region, who_df, sum, na.rm = TRUE)
case_lookup <- setNames(region_total_cases$cases_total, region_total_cases$region)

V(g_b)$total_cases <- case_lookup[V(g_b)$name]
V(g_b)$size <- rescale_to(log10(V(g_b)$total_cases + 1), to = c(15, 38))
V(g_b)$color <- "#2ca25f"
V(g_b)$label.cex <- 0.9
E(g_b)$width <- rescale_to(E(g_b)$similarity, to = c(1, 8))
E(g_b)$color <- "#636363"

png("figuras/rede_similaridade_regioes.png", width = 3000, height = 2200, res = 300)
plot(
	g_b,
	layout = layout_in_circle(g_b),
	main = "Estratégia B: Rede de similaridade entre regiões\nNós maiores = maior total acumulado de casos",
	vertex.label.color = "black",
	vertex.frame.color = NA,
	edge.curved = 0.15
)
legend(
	"topleft",
	legend = c("Região OMS", "Aresta de similaridade"),
	pt.bg = c("#2ca25f", NA),
	pch = c(21, NA),
	lty = c(NA, 1),
	col = c(NA, "#636363"),
	bty = "n",
	title = "Leitura de cores"
)
dev.off()

# -----------------------------------------------------------------------------
# Estratégia C: Rede de proximidade entre países (k-vizinhos mais próximos)
# -----------------------------------------------------------------------------

country_vars <- c("cases_per100k", "deaths_per100k")

country_df <- who_df[, c("country", "region", country_vars, "cases_total")]
country_df <- country_df[complete.cases(country_df[, country_vars]), ]
country_df <- country_df[order(country_df$cases_total, decreasing = TRUE), ]
country_df <- head(country_df, 120)

if (nrow(country_df) < 4) {
  stop("Dados insuficientes para a Estratégia C após limpeza.")
}

X <- scale(country_df[, country_vars])
d <- as.matrix(dist(X, method = "euclidean"))
diag(d) <- Inf

k <- min(3, nrow(country_df) - 1)

edge_list <- vector("list", nrow(country_df))
for (i in seq_len(nrow(country_df))) {
  nn <- order(d[i, ])[1:k]
  edge_list[[i]] <- data.frame(
    from = country_df$country[i],
    to = country_df$country[nn],
    dist = d[i, nn],
    stringsAsFactors = FALSE
  )
}

edges_c <- do.call(rbind, edge_list)

if (is.null(edges_c) || nrow(edges_c) == 0) {
  stop("Nao foi possivel gerar arestas para a Estratégia C.")
}

g_c <- graph_from_data_frame(edges_c, directed = FALSE, vertices = country_df)
g_c <- simplify(g_c, remove.multiple = TRUE, edge.attr.comb = "min")

com <- cluster_louvain(g_c)

comm_ids <- sort(unique(membership(com)))
palette_c <- setNames(rainbow(length(comm_ids)), as.character(comm_ids))
V(g_c)$color <- palette_c[as.character(membership(com))]
V(g_c)$size <- rescale_to(log10(V(g_c)$cases_total + 1), to = c(4, 18))
V(g_c)$label <- V(g_c)$name
V(g_c)$label.cex <- 0.3
E(g_c)$color <- "#bdbdbd"
E(g_c)$width <- 0.8

png("figuras/rede_proximidade_paises_knn.png", width = 3600, height = 2400, res = 300)
plot(
  g_c,
  layout = layout_with_fr(g_c),
  main = "Estratégia C: Rede k-NN entre paises (k=3)\nCores indicam comunidades por similaridade epidemiologica",
	vertex.label.color = "black",
  vertex.frame.color = NA,
  edge.curved = 0.05
)

legend_ids <- comm_ids
if (length(legend_ids) > 12) legend_ids <- legend_ids[1:12]

legend(
	"topleft",
	legend = paste("Comunidade", legend_ids),
	pt.bg = palette_c[as.character(legend_ids)],
	pch = 21,
	pt.cex = 1.1,
	cex = 0.7,
	bty = "n",
	title = "Cores dos nós"
)

if (length(comm_ids) > 12) {
	legend(
		"bottomleft",
		legend = "(Legenda truncada: primeiras 12 comunidades)",
		bty = "n",
		cex = 0.7,
		text.col = "#444444"
	)
}
dev.off()

# -----------------------------------------------------------------------------
# Estratégia D: Rede de perfil de letalidade e carga (pais-pais)
# -----------------------------------------------------------------------------

who_df$cfr <- ifelse(
	is.na(who_df$cases_total) | who_df$cases_total <= 0,
	NA_real_,
	who_df$deaths_total / who_df$cases_total
)

vars_d <- c("cases_per100k", "deaths_per100k", "cfr")

country_d <- who_df[, c("country", "region", "cases_total", vars_d)]
country_d <- country_d[complete.cases(country_d[, vars_d]), ]
country_d <- country_d[order(country_d$cases_total, decreasing = TRUE), ]
country_d <- head(country_d, 140)

if (nrow(country_d) >= 6) {
	Xd <- scale(country_d[, vars_d])
	dd <- as.matrix(dist(Xd, method = "euclidean"))
	diag(dd) <- Inf

	kd <- min(4, nrow(country_d) - 1)
	edge_list_d <- vector("list", nrow(country_d))

	for (i in seq_len(nrow(country_d))) {
		nn <- order(dd[i, ])[1:kd]
		edge_list_d[[i]] <- data.frame(
			from = country_d$country[i],
			to = country_d$country[nn],
			dist = dd[i, nn],
			stringsAsFactors = FALSE
		)
	}

	edges_d <- do.call(rbind, edge_list_d)

	if (!is.null(edges_d) && nrow(edges_d) > 0) {
		g_d <- graph_from_data_frame(edges_d, directed = FALSE, vertices = country_d)
		g_d <- simplify(g_d, remove.multiple = TRUE, edge.attr.comb = "min")

		reg_levels <- unique(V(g_d)$region)
		reg_palette <- setNames(rainbow(length(reg_levels)), reg_levels)

		V(g_d)$color <- reg_palette[V(g_d)$region]
		V(g_d)$size <- rescale_to(log10(V(g_d)$cases_total + 1), to = c(4, 18))
		V(g_d)$label <- V(g_d)$name
		V(g_d)$label.cex <- 0.4
		E(g_d)$width <- rescale_to(1 / (E(g_d)$dist + 1e-9), to = c(0.5, 2.4))
		E(g_d)$color <- "#bdbdbd"

		png("figuras/rede_perfil_letalidade_carga.png", width = 3600, height = 2400, res = 300)
		plot(
			g_d,
			layout = layout_with_fr(g_d),
			main = "Estratégia D: Rede de perfil de carga e letalidade\nArestas unem países com perfis semelhantes (casos/100k, obitos/100k, CFR)",
			vertex.label.color = "black",
			vertex.frame.color = NA,
			edge.curved = 0.03
		)
		legend(
			"topleft",
			legend = reg_levels,
			pt.bg = reg_palette[reg_levels],
			pch = 21,
			pt.cex = 1.1,
			cex = 0.75,
			bty = "n",
			title = "Regiões OMS"
		)
		dev.off()
	}
}

# -----------------------------------------------------------------------------
# Estratégia E: Rede bipartida Regiao-Indicador (assinatura regional)
# -----------------------------------------------------------------------------

indicator_vars <- c(
	"cases_total",
	"cases_per100k",
	"cases_7d",
	"cases_7d_per100k",
	"deaths_total",
	"deaths_per100k",
	"deaths_7d"
)

region_e <- aggregate(
	who_df[, indicator_vars],
	by = list(region = who_df$region),
	FUN = function(x) mean(x, na.rm = TRUE)
)

row.names(region_e) <- region_e$region
region_vals <- as.matrix(region_e[, -1])

z_mat <- scale(region_vals)
z_mat[is.na(z_mat)] <- 0

edge_e_idx <- which(z_mat >= 0.8, arr.ind = TRUE)

if (nrow(edge_e_idx) == 0) {
	# Fallback: para cada indicador, conecta a regiao com maior z-score
	edge_e_idx <- do.call(rbind, lapply(seq_len(ncol(z_mat)), function(j) {
		i <- which.max(z_mat[, j])
		c(i, j)
	}))
}

edges_e <- data.frame(
	from = row.names(z_mat)[edge_e_idx[, 1]],
	to = colnames(z_mat)[edge_e_idx[, 2]],
	weight = z_mat[edge_e_idx],
	stringsAsFactors = FALSE
)

if (nrow(edges_e) > 0) {
	g_e <- graph_from_data_frame(edges_e, directed = FALSE)

is_indicator <- V(g_e)$name %in% colnames(z_mat)
V(g_e)$color <- ifelse(is_indicator, "#7bccc4", "#ef6548")
V(g_e)$size <- ifelse(is_indicator, 12, 22)
V(g_e)$label.cex <- ifelse(is_indicator, 0.75, 0.9)
E(g_e)$width <- rescale_to(E(g_e)$weight, to = c(1, 6))
E(g_e)$color <- "#969696"

png("figuras/rede_bipartida_regiao_indicador.png", width = 3400, height = 2400, res = 300)
plot(
		g_e,
		layout = layout_with_fr(g_e),
		main = "Estratégia E: Rede bipartida Regiao-Indicador\nArestas representam indicadores com destaque regional (z-score >= 0.8)",
		vertex.label.color = "black",
		vertex.frame.color = NA,
		edge.curved = 0.07
	)
	legend(
		"topleft",
		legend = c("Indicador", "Região OMS"),
		pt.bg = c("#7bccc4", "#ef6548"),
		pch = 21,
		pt.cex = c(1.2, 1.6),
		bty = "n",
		title = "Cores dos nós"
	)
	dev.off()
}

###############################################################################
# 3) Resumo interpretativo para o relatório
###############################################################################

cat("\nResumo para discussão:\n")
cat("- Estratégia A mostra vínculo estrutural entre regiões e países com maior carga acumulada de casos.\n")
cat("- Estratégia B evidencia proximidade estatística entre regiões com perfis epidemiológicos semelhantes.\n")
cat("- Estratégia C revela agrupamentos de países (comunidades) por semelhança de indicadores per capita e recentes.\n")
cat("- Estratégia D compara países por perfil combinado de carga e letalidade (incluindo CFR).\n")
cat("- Estratégia E resume a assinatura epidemiológica das regiões via rede bipartida região-indicador.\n")
cat("\nArquivos gerados na pasta 'figuras'.\n")

