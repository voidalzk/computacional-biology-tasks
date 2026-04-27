library(igraph)
lines <- readLines('task-2/bio-yeast.mtx', warn = FALSE)
lines <- lines[!startsWith(lines, '%')][-1]
dados <- read.table(text = lines)
g <- simplify(graph_from_edgelist(as.matrix(dados[,1:2]), directed = FALSE))
comp <- components(g)
gc_idx <- which.max(comp[['csize']])
vertices_gc <- V(g)[comp[['membership']] == gc_idx]
gc <- induced_subgraph(g, vertices_gc)
deg <- degree(gc)
comm <- cluster_louvain(gc)
set.seed(42)
lay <- layout_with_fr(gc, niter = 1800)
labels <- rep(NA, vcount(gc))
hubs <- order(deg, decreasing = TRUE)[1:5]
labels[hubs] <- V(gc)top_hubs_grau_intermediacao_200ppi.png[hubs]
png('task-2/figuras/teste_gc_tmp.png', width=3200, height=2200, res=200, bg='white')
par(mar = c(1, 1, 3, 1))
plot(gc, layout=lay, vertex.size=2 + 8 * (deg - min(deg)) / (max(deg) - min(deg) + .Machine.eps), vertex.color=adjustcolor(hcl.colors(length(unique(membership(comm))), 'Dark 3')[membership(comm)], 0.92), vertex.label=labels, vertex.label.cex=0.7, vertex.label.color='black', vertex.frame.color=NA, edge.color=adjustcolor('gray35', 0.15), edge.width=0.45, main='bio-yeast (componente gigante)')
dev.off()
