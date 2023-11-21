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

```{math}
\def\bpi{\boldsymbol{\pi}}
\def\bzero{\boldsymbol{0}}
\def\be{\boldsymbol{e}}
```

# Non-equilibrium steady states

The equilibrium formulation assumes detailed balance. However, when the state-transition diagram of a receptor model includes cycles, the steady-state probability distribution may not satisfy detailed balance.

Let us consider again the following state-transition diagram. 

```{code-cell}
var('a12, a21, a13, a31, a23, a32, a34, a43')
d = {1: {2:a12, 3:a13}, 2: {1:a21, 3:a23}, 3: {2:a32, 1:a31, 4:a34}, 4: {3:a43}};
G = DiGraph(d,weighted=True)
vertex_positions = {1: (0, 0), 2: (1, 1.41), 3: (2, 0), 4: (4,0)}
G.plot(figsize=8,edge_labels=True,pos=vertex_positions,graph_border=True)
```

## Generator matrix and Laplacian 

The generator matrix `Q` for the Markov chain associated to `G` can be constructed from the weighted adjacency matrix `A`, given by 
```{code-cell}
A = G.weighted_adjacency_matrix()
print(A)
```

```{code-cell}
Q = A - diagonal_matrix(sum(A.T))
print(Q)
```

Equivalently, the generator matrix `Q` is the opposite of the Laplacian matrix `L` (
`Q = -G.laplacian_matrix()` (not shown).  `Q` is referred to as the generator matrix for the Markov chain because the probability distribution {math}`\bpi^T` solves {math}`d\bpi^T/dt = \bpi^T Q`.  

```{note}
The probability distribution {math}`\bpi^T` is a row vector.
```

The steady-state probability distribution solves {math}`\bpi^T Q = \bzero` subject to {math}`\bpi^T \be = 1` (conservation of probability).  In this expression, {math}`\be` is a commensurate column vector of 1s. 

The following code confirms that each row of `Q` sums to zero.
```{code-cell}
e = matrix([1,1,1,1]).T
print(Q*e)
```
To see that the columns of `Q` do not sum to zero, we evaluate {math}`\be^T Q`.

```{code-cell}
print(e.T*Q)
```

## Symbolic solution 

## Condition on cycles 
