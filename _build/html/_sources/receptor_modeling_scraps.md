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


# Non-equilibrium steady states

The equilibrium formulation assumes detailed balance. However, when the state-transition diagram of a receptor model includes cycles, the steady-state probability distribution may not satisfy detailed balance.

Let us consider again the following state-transition diagram. 

```{code-cell}
var('a12, a21, a13, a31, a23, a32, a34, a43')
d = {1: {2:a12, 3:a13}, 2: {1:a21, 3:a23}, 3: {2:a32, 1:a31, 4:a34}, 4: {3:a43}};
G = DiGraph(d,weighted=True)
G.plot(figsize=8,edge_labels=True,pos=vertex_positions,graph_border=True)
```

The weighted adjacency matrix for `G` is 

```{code-cell}
A = G.weighted_adjacency_matrix()
A
```

## Generator matrix and Laplacian 

The generator matrix `Q` for the Markov chain associated to `G` can be constructed from the weighted adjacency matrix `A` as folows.


```{code-cell}
Q = A - diagonal_matrix(sum(A.T))
```

Equivalently, the generator matrix `Q` is the opposite of the Laplacian matrix `L` ({math}`Q=-L`). The Laplacian of `G` is
```{code-cell}
L = G.laplacian_matrix()
L
```






## Scraps

The following code defines `e` to be column vector of ones.  This is used to show that each row of `Q` sums to zero.
```{code-cell}
e = matrix([1,1,1,1]).T
```

```{code-cell}
print(Q*e)
```

Multiplying on the left by the transpose of `e`, given by `e.T`, we see that each column of `Q` does not sum to zero.

```{code-cell}
print(e.T*Q)
```

