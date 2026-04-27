import igraph as ig
import matplotlib.pyplot as plt
import os

arquivo_dataset = "bio-yeast.mtx"
arquivo_limpo = "bio-yeast-limpo.edges"

with open(arquivo_dataset, 'r') as f_in, open(arquivo_limpo, 'w') as f_out:
    for linha in f_in:
        if not linha.startswith('%'):
            f_out.write(linha)

grafo = ig.Graph.Read_Ncol(arquivo_limpo, directed=False)
grafo.simplify(multiple=True, loops=True)

print("Direcionado:", grafo.is_directed())
print("Ordem (Vértices):", grafo.vcount())
print("Tamanho (Arestas):", grafo.ecount())
print("Densidade:", grafo.density())
print("Diâmetro:", grafo.diameter())
print("Transitividade (Clustering):", grafo.transitivity_undirected())

graus = grafo.degree()
intermediacoes = grafo.betweenness()

top_5_grau = sorted(range(len(graus)), key=lambda i: graus[i], reverse=True)[:5]
top_5_intermediacao = sorted(range(len(intermediacoes)), key=lambda i: intermediacoes[i], reverse=True)[:5]

print("Top 5 nós (Grau):", [grafo.vs[i]["name"] for i in top_5_grau])
print("Top 5 nós (Intermediação):", [grafo.vs[i]["name"] for i in top_5_intermediacao])

comunidades = grafo.community_multilevel()

fig, ax = plt.subplots(figsize=(12, 12))
layout = grafo.layout_fruchterman_reingold()

tamanhos_nos = [max(5, grau * 0.5) for grau in graus]

ig.plot(
    comunidades,
    palette=ig.ClusterColoringPalette(len(comunidades)),
    vertex_size=tamanhos_nos,
    layout=layout,
    target=ax
)

plt.show()

if os.path.exists(arquivo_limpo):
    os.remove(arquivo_limpo)