ANÁLISE DE REDES COMPLEXAS COM IGRAPH
O objetivo deste exercício é aplicar conceitos de Teoria dos Grafos para extrair insights a partir de uma rede real, utilizando programação tanto para manipulação de dados quanto para análise estrutural. A atividade pode ser realizada individualmente ou em grupos de até três alunos.

Cada grupo deverá selecionar e baixar um dataset de rede disponível no repositório Network Repository (https://networkrepository.com/). Este repositório reúne redes de diferentes domínios disponibilizadas pela comunidade acadêmica.

Recomendação: priorizar redes das categorias Social Network ou Biology.

Restrição: selecionar redes com 200 a 2000 nós para garantir que o processamento e a visualização sejam viáveis em seus computadores.

Formatos comuns: .edges (lista de arestas) ou .mtx (Matrix Market).

# O DESAFIO: DESENVOLVER UM SCRIPT DE ANÁLISE DE REDES

Utilizando a biblioteca igraph (em R ou Python), desenvolvam um script que execute as seguintes etapas:

Etapa 1. Importação e Limpeza

– Importar o dataset e convertê-lo em um objeto de grafo.

– Identificar se a rede é direcionada ou não-direcionada.

– Realizar a simplificação da rede: Remover loops (arestas que ligam um nó a ele mesmo) e arestas múltiplas.

Etapa 2. Caracterização Topológica

Calculem e apresentem os seguintes indicadores globais da rede:

– Ordem e Tamanho (número total de vértices e arestas)

– Densidade (proporção entre conexões existentes e possíveis)

– Diâmetro (maior distância geodésica entre dois nós)

– Transitividade (coeficiente de agrupamento da rede; clustering coefficient)

Etapa 3. Identificação de Atores-chave (Hubs)

Calculem métricas de centralidade para todos os nós e identifique os 5 nós mais importantes segundo o grau (degree) e intermediação (Betweenness). Discutir o papel estrutural desses nós (por exemplo, conectores locais vs. pontes globais).

Etapa 4. Detecção de Comunidades e Visualização

Apliquem um algoritmo de detecção de comunidades (ex: Louvain ou Walktrap) e gerem uma visualização onde o tamanho do nó seja proporcional ao seu Grau, a cor do nó represente a comunidade à qual ele pertence, e o layout seja de Fruchterman-Reingold ou equivalente.

# RELATÓRIO DE ENTREGA

Além do código comentado, entreguem uma breve análise (1 a 2 páginas) indicando o tipo da rede escolhida, descrevendo os Hubs e as comunidades encontradas. Para identificação do tipo de rede, considere métricas como densidade e diâmetro, as quais podem auxiliar na classificação da rede (e.g. “mundo pequeno” ou outro tipo de rede). Para análise dos Hubs, indique se os nós com maior grau são os mesmos com maior Betweenness. Se houver diferença, explique o que isso revela sobre a estrutura da rede escolhida. Para interpretação das comunidades, observe a representação da rede em um grafo. Olhando para a imagem, observe se as comunidades encontradas fazem sentido visual. O que elas poderiam representar no contexto do seu dataset (por exemplo, se for uma rede social, seriam grupos de amigos? Se for biológica, seriam genes com funções similares?).

Dica: Em R: Use read.table() ou Matrix::readMM() para ler os dados e graph_from_data_frame() para criar o grafo.