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
for n in range(2,7):
  G=graphs.PathGraph(n)
  G.show(figsize=4,graph_border=True,title=f'P{n}')
```


