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

# Example graphs

(example_graphs:path_graph)=
## The path graph {math}`P_n`

A `path graph` is a graph whose {math}`n` vertices can be listed in the order {math}`0, 1, \ldots , n-1` such that the edges are {math}`(i, i+1)` for {math}`i = 0, 1 \ldots , n-2`.  For example, the undirected graphs {math}`P_2` through {math}`P_6` are shown below.

```{code-cell}
:tags: [hide-input]
for n in range(2,7):
  G=graphs.PathGraph(n)
  G.show(figsize=4,graph_border=True,title=f'P{n}')
```

(example_graphs:cycle_graph)=
## The cycle graph $C_n$

The  `cycle graph` $C_n$ is a graph whose $n$ vertices can be listed in the order $0, 1, \ldots , n$ such that the edges are $(i, i+1)$ for $i = 0, 1, \ldots , n-2$, 
and also $(0,n-1)$.  For example, $C_3$, $C_5$, and $C_{12}$ are 

```{code-cell}
:tags: [hide-input]

for n in [3,5,12]:
    graphs.CycleGraph(n).show(figsize=3,title='K%s:' %(n))
```

(example_graphs:complete_graph)=
## The complete graph $K_n$ 

The   `complete graph` $K_n$ is a graph with $n$ vertices and $|E(K_n)|=\binom{n}{2}$ edges.  
For example, $K_2$, $K_4$ and $K_7$ are

```{code-cell}
:tags: [hide-input]

for n in [2,4,7]:
    graphs.CompleteGraph(n).show(figsize=3,title='C%s:' %(n))
```

(example_graphs:hypercube_graph)=
## The hypercube graph $Q_n$

The `hypercube graph`  $Q_n$ has $|V(Q_n)|=2^n$ vertices and $|E(Q_n)|=2^{n-1}n $ edges.  For example, $Q_2$, $Q_3$ and $Q_5$ are 

```{code-cell}
:tags: [hide-input]

for n in [2,3,5]:
    if n>3:
        graphs.CubeGraph(n).show(figsize=3,vertex_size=20,vertex_labels=false,title='Q%s:' %(n))
    else:
        graphs.CubeGraph(n).show(figsize=3,title='Q%s:' %(n))
```

Here is a 3D image of $Q_4$:
```{code-cell}
:tags: [hide-input]

n=4
graphs.CubeGraph(n).show3d(title='Q%s:' %(n))
```
