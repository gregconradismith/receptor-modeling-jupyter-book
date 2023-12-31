---
jupytext:
  cell_metadata_filter: -all
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.10.3
kernelspec:
  display_name: SageMath 10
  language: sage
  name: sage-10.0
---

# States and transitions

The process of receptor modeling often begins by specifying the molecular conformations (states) to be considered and the transitions between these states.  

## Undirected graphs as state-transition diagrams

For example, the following graph may represent a receptor model that has four states:

```{code-cell}
G_undirected = Graph({1: [2, 3], 2: [3], 3: [4]})
vertex_positions = {1: (0, 0), 2: (1, 1.41), 3: (2, 0), 4: (4,0)}
G_undirected.plot(figsize=8,pos=vertex_positions,graph_border=True)
```

```{code-cell}
G_undirected = Graph({1: [2, 3], 2: [3], 3: [4]})
vertex_positions = {1: (0, 0), 2: (1, 1.41), 3: (2, 0), 4: (4,0)}
G_undirected.plot(figsize=8,pos=vertex_positions,graph_border=True,talk=True)
```

The graph `G_undirected` is constructed by calling the Sagemath command [`Graph()`](https://doc.sagemath.org/html/en/reference/graphs/sage/graphs/graph.html#supported-formats) with a dictionary that associates neighbors to each vertex.  The vertices of the graph `G` are the integers 1, 2, 3, and 4.  The method `plot()` shows the graph `G_undirected` using a dictionary `vertex_positions` that specifies the locations of each vertex.

The adjacency matrix of this graph is
```{code-cell}
G_undirected.adjacency_matrix()
```

The incidence matrix of this graph is
```{code-cell}
G_undirected.incidence_matrix()
```

```{note}
The graphs used here to represent receptor states and transitions will be both connected and simple.
* _Connected_: at least one path joins every pair of vertices.
* _Simple_: no loops or multiple edges.
```

## Directed graphs as state-transition diagrams 

In the context of receptor modeling, the undirected graph above is interpreted as a short-hand for the following _directed_ graph.

```{code-cell}
G_directed = G_undirected.to_directed()
G_directed.plot(figsize=8,edge_labels=True,pos=vertex_positions,graph_border=True)
```

```{code-cell}
G_directed = G_undirected.to_directed()
G_directed.plot(figsize=8,edge_labels=True,pos=vertex_positions,graph_border=True,talk=True)
```

The method `to_directed()` produces `G_directed` as the _symmetric_ digraph associated to `G_undirected`, in which adjacent vertices are  connected in both directions.

```{note}
Receptor state-transition diagrams  will always be symmetric directed graphs, that is, for every edge from vertex `i` to vertex `j`, there is also an edge from vertex `j` to vertex `i`.  Thus, the state-transtion diagrams for a receptor model may, for simplicity, be illustrated as an undirected graph.
```

## Transition rate constants 

In the context of receptor modeling, state-transition diagrams are usually _weighted_, as shown here.

```{code-cell}
var('a12, a21, a13, a31, a23, a32, a34, a43')
d = {1: {2:a12, 3:a13}, 2: {1:a21, 3:a23}, 3: {2:a32, 1:a31, 4:a34}, 4: {3:a43}};
G = DiGraph(d,weighted=True)
G.plot(figsize=8,edge_labels=True,pos=vertex_positions,graph_border=True)
```

In the code above, a directed graph `G` is constructed by calling the Sagemath command [`DiGraph()`](https://doc.sagemath.org/html/en/reference/graphs/sage/graphs/digraph.html#methods).  The input argument `d` is a [Python dictionary](https://doc.sagemath.org/html/en/thematic_tutorials/tutorial-programming-python.html) that assigns out-neighbors to each vertex and corresponding edge labels.
The edge labels are not _strings_, but _symbolic variables_ defined using Sagemath's `var` command.
For example, the symbolic variable `a12` stands for the rate of transition between state 1 and 2. 

Because these rate constants are symbolic variables, Sagemath will evaluate expressions such as
```{code-cell}
f = a12 * (a21 + a13)^2 / a13
f.expand()
```

One reason for using symbolic variables is that we can produce symbolic expressions important quantities using [module for graphs and digraphs](https://doc.sagemath.org/html/en/reference/graphs/index.html) available in `Sagemath`. For example, the weighted adjacency matrix for `G` is 

```{code-cell}
G.weighted_adjacency_matrix()
```

The Laplacian of `G` is
```{code-cell}
G.laplacian_matrix()
```
This matrix is sometimes referred to as the `combinatorial Laplacian matrix` of the weighted directed graph `G`.



