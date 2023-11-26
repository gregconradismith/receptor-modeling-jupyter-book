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

A `path graph` is a graph whose {math}`n` vertices can be listed in the order {math}`0, 1, \ldots , n-1` such that the edges are {math}`(i, i+1)` for {math}`i = 0, 1 \ldots , n-2`.  For example, the undirected graphs {math}`P_3` and {math}`P_6` are shown below.

```{code-cell}
:tags: [hide-input]
for n in [3,6]:
  G=graphs.PathGraph(n)
  G.show(figsize=4,graph_border=True,title=f'P{n}')
```

(example_graphs:cycle_graph)=
## The cycle graph {math}`C_n`

The  `cycle graph` {math}`C_n` is a graph whose $n$ vertices can be listed in the order {math}`0, 1, \ldots , n` such that the edges are {math}`(i, i+1)` for {math}`i = 0, 1, \ldots , n-2`, 
and also {math}`(0,n-1)`.  For example, {math}`C_3`, {math}`C_5`, and {math}`C_{12}` are 

```{code-cell}
:tags: [hide-input]

for n in [3,5,12]:
    graphs.CycleGraph(n).show(figsize=3,title='C%s:' %(n))
```

(example_graphs:complete_graph)=
## The complete graph {math}`K_n` 

The   `complete graph` {math}`K_n` is a graph with {math}`n` vertices and {math}`|E(K_n)|=\binom{n}{2}` edges.  
For example, {math}`K_2`, {math}`K_4`, and {math}`K_7` are

```{code-cell}
:tags: [hide-input]

for n in [2,4,7]:
    graphs.CompleteGraph(n).show(figsize=3,title='K%s:' %(n))
```

(example_graphs:hypercube_graph)=
## The hypercube graph {math}`Q_n`

The `hypercube graph`  {math}`Q_n` has {math}`|V(Q_n)|=2^n` vertices and {math}`|E(Q_n)|=2^{n-1}n` edges.  For example, {math}`Q_2` and {math}`Q_5` are 

```{code-cell}
:tags: [hide-input]

for n in [2,5,9]:
    if n>3:
        graphs.CubeGraph(n).show(figsize=3,vertex_size=20,vertex_labels=false,title='Q%s:' %(n))
    else:
        graphs.CubeGraph(n).show(figsize=3,title='Q%s:' %(n))
```

Here is a 3D image of {math}`Q_4`:
```{code-cell}
:tags: [hide-input]

n=4
graphs.CubeGraph(n).show3d(title='Q%s:' %(n))
```

Here is a 3D image of {math}`Q_3` and {math}`Q_4`:
```{code-cell}
:tags: [hide-input]

for n in [3,4]:
    graphs.CubeGraph(n).show3d(title='Q%s:' %(n))
```
