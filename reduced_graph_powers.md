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

# Reduced graph powers

The Sagemath `Graph()` object will be used for the topology of receptor models.

## Simple graphs

For a simple graph {math}`G=(V,E)`, we will sometimes write {math}`|V(G)|` and {math}`|E(G)|`
for the number of vertices and edges of {math}`G`.  We assume the graph {math}`G` is connected, so the Betti number is given by

```{math}
\beta(G)=|E(G)|-|V(G)|+1 \, .
```

When the graph {math}`G` is understood from context, we will often use {math}`v(G)` and {math}`e(G)`, or simply {math}`v` and {math}`e`, for the number of vertices and edges of {math}`G`.

The number of vertices in the reduced graph power {math}`G^{(k)}` is

```{math}
|V(G^{(k)})|=\binom{v+k-1}{k}
```

where {math}`v=|V(G)|`. The number of edges in the reduced graph power {math}`G^{(k)}` is

```{math}
:label: number_of_edges
|E(G^{(k)})|=e\binom{v+k-2}{k-1}
```

where {math}`e=|E(G)|`.
When the graph {math}`G` is understood, will sometimes abbreviate the number of vertices and edges in its reduced graph power by {math}`v^{(k)} = |V(G^{(k)})|` and {math}`e^{(k)} = |E(G^{(k)})|`.

## Common named graphs

### The path graph {math}`P_n`

A `path graph` is a graph whose {math}`n` vertices can be listed in the order {math}`0, 1, \ldots , n-1` such that the edges are {math}`(i, i+1)` for {math}`i = 0, 1 \ldots , n-2`.  For example,

```{code-cell}
for n in [2,4,7]:
    graphs.PathGraph(n).show(figsize=3,graph_border=True,title='P%s:' %(n))
```

Evidently, {math}`|V(P_n)|=n`, {math}`|E(P_n)|=n-1`, and {math}`\beta(P_n)=0`. The number of vertices of the reduced graph power {math}`P_n^{(k)}` is {math}`|V(P_n^{(k)})|=\binom{n+k-1}{k}`. The number of edges is {math}`|E(P_n^{(k)})|=(n-1)\binom{n+k-2}{k-1}`. Consequently, the Betti number of the path graph on {math}`n` vertices is

```{math}
\begin{equation}
\beta(P_n^{(k)}) = (n-1)\binom{n+k-2}{k-1}-\binom{n+k-1}{k}+1 \, .
\end{equation}
```

Specializing to the case of a dimer ({math}`k=2`), we have {math}`|V(P_n^{(2)})|=\binom{n+1}{2}`,  {math}`|E(P_n^{(2)})|=n(n-1)`, and {math}`\beta(P_n^{(2)}) = \binom{n-1}{2}`.

```{math}
:label: beta_for_path_graph
\beta(P_n^{(k)}) = (n-1)\binom{n+k-2}{k-1}-\binom{n+k-1}{k}+1 \, .
```

### The cycle graph {math}`C_n`

The  `cycle graph` {math}`C_n` is a graph whose {math}`n` vertices can be listed in the order {math}`0, 1, \ldots , n` such that the edges are {math}`(i, i+1)` for {math}`i = 0, 1, \ldots , n-2`, and also {math}`(0,n-1)`.  For example, {math}`C_5` is

```{code-cell}
:tags: [hide-input]
for n in [3,5,12]:
    graphs.CycleGraph(n).show(figsize=3,title='K%s:' %(n))
```

For the cycle graph {math}`|E(C_n)|=n` and {math}`\beta(C_n)=1`.

The number of vertices of the reduced graph power {math}`C_n^{(k)}` is {math}`|V(C_n^{(k)})|=\binom{n+k-1}{k}`. The number of edges is {math}`|E(C_n^{(k)})|=n\binom{n+k-2}{k-1}`. The Betti number is

```{math}
\beta(C_n^{(k)}) = n\binom{n+k-2}{k-1}-\binom{n+k-1}{k}+1 \, .
```

Specializing to the case of a dimer ({math}`k=2`), we have {math}`v^{(2)}=\binom{n+1}{2}`,  {math}`e^{(2)}=n^2`, and

```{math}
\beta(C_n^{(2)}) =  \binom{n}{2}+1 
= n(n-1)/2+1 \,.
```

### The complete graph {math}`K_n`

The   `complete graph` {math}`K_n` is a graph with {math}`n` vertices and {math}`|E(K_n)|=\binom{n}{2}` edges.
For example, $K_5$ is

```{code-cell}
:tags: [hide-input]

for n in [2,4,7]:
    graphs.CompleteGraph(n).show(figsize=3,title='C%s:' %(n))
```

For the complete graph, {math}`\beta(K_n)=\binom{n-1}{2}=(n-1)(n-2)/2`.

The number of vertices of the reduced graph power ${math}`K_n^{(k)}` is {math}`|V(K_n^{(k)})|=\binom{n+k-1}{k}`, the number of edges is {math}`|E(K_n^{(k)})|=\binom{n}{2}\binom{n+k-2}{k-1}`, and the Betti number is

```{math}
\beta(K_n^{(k)}) = \binom{n}{2}\binom{n+k-2}{k-1}-\binom{n+k-1}{k}+1 \, .
```

For a dimer ({math}`k=2`), we have {math}`v^{(2)}=\binom{n+1}{2}`,  {math}`e^{(2)}=\binom{n}{2} n = n^2(n-1)/2`, and

```{math}
\beta(K_n^{(2)}) = (n+1)(n-1)(n-2)/2 \, .
```

### The hypercube graph {math}`Q_n`

The `hypercube graph`  {math}`Q_n` has {math}`|V(Q_n)|=2^n` vertices and {math}`|E(Q_n)|=2^{n-1}n` edges.  For example, {math}`Q_4` is

```{code-cell}
:tags: [hide-input]

for n in [2,4,5]:
    if n>3:
        graphs.CubeGraph(n).show(figsize=3,vertex_size=20,vertex_labels=false,title='Q%s:' %(n))
    else:
        graphs.CubeGraph(n).show(figsize=3,title='Q%s:' %(n))
```

```{code-cell}
:tags: [hide-input]

n=3
graphs.CubeGraph(n).show3d(title='Q%s:' %(n))
```

For the hypercube graph, {math}`\beta(Q_n)= 2^{n-1} (n-2)+1`.

The number of vertices of the reduced graph power {math}`Q_n^{(k)}` is {math}`|V(Q_n^{(k)})|=\binom{2^n+k-1}{k}`.  The number of edges is {math}`|E(Q_n^{(k)})|=2^{n-1} n \binom{2^n+k-2}{k-1}`. The Betti number is

```{math}
\beta(Q_n^{(k)}) =2^{n-1} n \binom{2^n+k-2}{k-1}-\binom{2^n+k-1}{k}+1 \, .
```

The reduced graph product of two hypercube graphs, has {math}`|V(Q_n^{(2)})|=2^{n-1} (2^n+1)` vertices,  {math}`|E(Q_n^{(2)})| = 2^{2n-1} n` edges, and Betti number

```{math}
\beta(Q_n^{(2)}) =  2^{2n-1} (n-1) - 2^{n-1}+1 \, . 
```

In particular, {math}`\beta(Q_2^{(2)}) = 7` and  {math}`\beta(Q_3^{(2)}) = 61`.

An excellent resource on hypergraphs is {cite}`Harary1988`
[{PDF}](https://www.sciencedirect.com/science/article/pii/0898122188902131).


## References 

```{bibliography}
:filter: docname in docnames
```

