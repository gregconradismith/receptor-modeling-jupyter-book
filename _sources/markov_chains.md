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
# Markov chain receptor models 

## Symbolic vertex and edge labels

Recall the [three state model](receptor_modeling_ligands:three_state_model) discussed above ([Ligands](receptor_modeling_ligands)).  Using generic notation for the [path graph](example_graphs:path_graph) on three vertices, 

```{code-cell}
var('p0 p1 p2 a01 a10 a12 a21')
G=graphs.PathGraph(3).to_directed()
G.relabel({0:p0,1:p1,2:p2})
G.set_edge_label(p0,p1,a01)
G.set_edge_label(p1,p0,a10)
G.set_edge_label(p1,p2,a12)
G.set_edge_label(p2,p1,a21)
G.show(figsize=4,edge_labels=True,talk=True)
```

We assign symbolic variables to the vertices and edges of $G$ because this allows us to produce symbolic expressions important quantities using methods available in `Sagemath`'s [module for graphs and digraphs](https://doc.sagemath.org/html/en/reference/graphs/index.html). For example, the weighted adjacency matrix associated with graph $G$ above is

```{code-cell}
A = G.weighted_adjacency_matrix()
```

The combinatorial Laplacian matrix is $L=D-A$ where $D$ is a diagonal matrix given by the column sum of $A$.  

```{code-cell}
L = diagonal_matrix(sum(A.T))-A
show(L)
```

```{note}
In the code block above, I prefer to write ```sum(A.T)``` rather than ```sum(A.columns())```, but these are equal.
```

```{code-cell}
sum(A.T) == sum(A.columns())
```

## Generator matrix of Markov chain 

The Markov chain with states and transistions of $G$ has a generator matrix $Q$ that is the opposite (additive inverse) of the Laplacian matrix ($L$). That is, $Q=-L=A-D$. We define a function that calculates the symbolic generator matrix $Q$ from the weighted adjacency matrix $A$.

```{code-cell}
def generator(A):
    return A-diagonal_matrix(sum(A.T))
```

Calling the function ```generator()``` gives the result we expect:

```{code-cell}
Q = generator(A)
show(Q)
```

Note that the sums of the columns are zero, reflecting the conservation of probability.

```{code-cell}
sum(Q.columns()) == sum(Q.T)
```

## Markov chain tree theorem (Hill's diagrammatic method)

Using the Markov chain tree theorem \(CHECK\), it is straightforward to find symbolic expressions for the steady\-state probabilities of each state \(p0,p1,p2\).  The relative probabilities are:

```{code-cell}
z0 = Q[[1,2],[1,2]].determinant().simplify_full()
z1 = Q[[0,2],[0,2]].determinant().simplify_full()
z2 = Q[[0,1],[0,1]].determinant().simplify_full()
print(f'[ {z0} : {z1} : {z2}]')
```

The normalized probabilities are:

```{code-cell}
zT = z0+z1+z2
p0 = z0/zT
p1 = z1/zT
p2 = z2/zT
show(table([[f'{p0=}'],[f'{p1=}'],[f'{p2=}']]))
```

This calculation can be peformed for a Markov chain whose topology is a simple directed graph (connected, no loops) that is also symmetric.

```{code-cell}
def hill_diagramatic_method(Q):
    n = Q.nrows()
    if Q.ncols() != n:
        raise ValueError
    z = [0]*n
    for i in range(n):
        a = [ j for j in range(n) ]
        a.remove(i)
        z[i] = Q[a,a].determinant().simplify_full()
    return z

G=graphs.HouseGraph().to_directed()
# put values on each edge ... 
zz = hill_diagramatic_method(Q)
print(zz)
```

```{code-cell}
# define this? G.weighted_laplacian_matrix()
```

