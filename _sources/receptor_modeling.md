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


# Receptor models


Quantitative pharmacologists construct models of relationship between ligand concentration and the fraction of cell surface receptors in each of several molecular conformations. These models give insight into the action of natural ligands and drugs on receptor-mediated cell responses.  A good starting point for review of receptor modeling is 

* Kenakin, T. P. (2018). [A Pharmacology Primer: Techniques for More Effective and Strategic Drug Discovery.](https://www.amazon.com/Pharmacology-Primer-Techniques-Effective-Strategic/dp/0128139579/) Academic Press. {cite}`Kenakin2018`

An introction to this topic for undergraduate neuroscienc majors is __Ligands, Receptors, and Rate Laws__, Chapter 4 of [Cellular Biophysics \& Modeling: a primer on the computational biology of excitable cells](https://www.cambridge.org/core/books/cellular-biophysics-and-modeling/0C728F4C44D89D8F3BA62E41A0D7336F). Cambridge University Press, 2019.

## States and transitions

The process of receptor modeling often begins by specifying the molecular conformations (states) to be considered and the transitions between these states.  For example, the following graph may represent a receptor model that has four states:

```{code-cell}
G_undirected = Graph({1: [2, 3], 2: [3], 3: [4]})
vertex_positions = {1: (0, 0), 2: (1, 1.41), 3: (2, 0), 4: (4,0)}
G_undirected.plot(figsize=8,pos=vertex_positions,graph_border=True)
```

This graph `G_undirected` is constructed by calling the Sagemath command [`Graph()`](https://doc.sagemath.org/html/en/reference/graphs/sage/graphs/graph.html#supported-formats) with a dictionary that associates neighbors to each vertex.  The vertices of the graph `G` are the integers 1, 2, 3, and 4.  The method `plot()` shows the graph `G_undirected` using a dictionary `vertex_positions` that specifies the locations of each vertex.

The adjacency matrix of this graph is
```{code-cell}
A = G_undirected.adjacency_matrix()
```

The incidence matrix of this graph is
```{code-cell}
E = G_undirected.incidence_matrix()
```

```{note}
The graphs used here to represent receptor states and transitions will be both connected and simple.
* _Connected_: at least one path joins every pair of vertices.
* _Simple_: no loops or multiple edges.
```

In the context of receptor modeling, the undirected graph above is interpreted as a short-hand for the following _directed_ graph.

```{code-cell}
G_directed = G_undirected.to_directed()
G_directed.plot(figsize=8,edge_labels=True,pos=vertex_positions,graph_border=True)
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

The weighted adjacency matrix for `G` is 

```{code-cell}
G.weighted_adjacency_matrix()
```

The Laplacian of `G` is:
```{code-cell}
G.laplacian_matrix()
```
This matrix is sometimes referred to as the `combinatorial Laplacian matrix` of the weighted directed graph `G`.




The generator matrix `Q` for the Markov chain associated to `G` can be constructed from the weighted adjacency matrix `A` as folows.
```{code-cell}
Q = A - diagonal_matrix(sum(A.T))
```

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




(receptors:three_state_model)=
## A three-state receptor model

As a simple example, consider a receptor model with three states arranged as follows.

```{code-cell}
G = DiGraph({0: {1:'a01'}, 1: {0:'a10', 2:'a12'}, 2: {1:'a21'}})
pos = {0: (0, 0), 1: (1, 0), 2: (2, 0)}
G.plot(figsize=4,edge_labels=True,pos=pos,graph_border=True)
```

```{code-cell}
G = DiGraph({'R': {'RL':'kap*L'}, 'RL': {'R':'kam', 'RLL':'kbp*L'}, 'RLL': {'RL':'kbm'}})
pos = {'R': (0, 0), 'RL': (1, 0), 'RLL': (2, 0)}
G.plot(figsize=4,edge_labels=True,pos=pos,graph_border=True,vertex_size=1000)
```


    
    
When both forward and reverse transitions are explicit, the state-transition diagram has the topology of a symmetric directed version of the [path graph](example_graphs:path_graph) with 3 vertices.

Here is the (undirected) path graph {math}`P_3`:

```{code-cell}
G = Graph({0: [1], 1: [2]})
G.plot(figsize=4)
```

The symmetric directed version is

```{code-cell}
:tags: ["hide-input"]
G = DiGraph({0: [1], 1: [0,2], 2: [1]})
G.plot(figsize=4)
```

Our practice will be to define symbolic variables and put these on the vertices and edges.

```{code-cell}
var('a b c kappa_b_plus kappa_b_minus kappa_c_plus kappa_c_minus')
G = DiGraph([[a,b,c],[(a,b),(b,a),(b,c),(c,b)]])
G.set_edge_label(a,b,kappa_b_plus)
G.set_edge_label(b,a,kappa_b_minus)
G.set_edge_label(b,c,kappa_c_plus)
G.set_edge_label(c,b,kappa_c_minus)
G.plot(figsize=4,edge_labels=True)
```

In the state-transition diagram shown above, {math}`\kappab` and {math}`\kappac` are dimensionless equilibrium constants, {math}`\kappabstar` and {math}`\kappacstar` are association constants with physical dimension of inverse concentration, and {math}`x` is ligand concentration.  The solid harpoons indicate the forward reaction direction.  For example, the reaction labelled {math}`\kappab` has {math}`a` as reactant and {math}`b` as product; consequently, increasing {math}`\kappab` decreases the equilibrium probability (relative fraction) of state {math}`a` and increases the probability of state {math}`b`.
The three states of Equation XXX are labelled so that the reactant comes before the product in dictionary order ({math}`a` to {math}`b` to {math}`c`).  The subscript of the equilibrium constants {math}`\kappab` and {math}`\kappac` are chosen to match the label of the reaction products.


For an isolated monomer with a state-transition diagram given by Equation XXX, the probability of state {math}`i` is given by
{math}`\pi_i = z_i / z_T` where  {math}`z_T= \textstyle \sum_i z_i`,
{math}`z_a = 1`,
{math}`z_b = \kappab = \kappabstar x`, and
{math}`z_c =\kappab \kappac = \kappabstar \kappacstar x^2`. That is,

\begin{equation}
\pi_a =  \frac{1}{1+ \kappabstar x  +  \kappabstar \kappacstar x^2} \, ,  \quad \pi_b =  \frac{\kappabstar x}{1+ \kappabstar x +  \kappabstar  \kappacstar x^2}   \quad \mbox{and}  \quad \pi_c = \frac{\kappabstar  \kappacstar x^2 }{1+ \kappabstar x +  \kappabstar \kappacstar x^2 }  \, .
\end{equation}

It is helpful to  present this set of rational functions using the following compact notation:
\begin{equation}
 [ \pi_a  :  \pi_b :  \pi_c ] = [1 :\kappab : \kappab \kappac ]  = [1 : \kappabstar x :\kappabstar \kappacstar x^2  ] \,  .
\end{equation}
In expressions of this kind,  it is understood that
\begin{equation}
[ x_1 \! : \! x_2 : \! \cdots \! : \! x_n ] = [ \lambda x_1 \! : \! \lambda  x_2 : \! \cdots \! : \! \lambda  x_n ]
\end{equation}
for any {math}`\lambda \neq 0`. Furthermore, {math}`\lambda = 1/\sum_i x_n` gives the probability distribution {math}`\pi = (\pi_1, \pi_2, \ldots, \pi_n)` where {math}`1=\sum_i \pi_i`.

One reason for using symbolic variables is that we can produce symbolic expressions important quantities using [module for graphs and digraphs](https://doc.sagemath.org/html/en/reference/graphs/index.html) available in `Sagemath`. For example, the weighted adjacency matrix associated with graph {math}`G$ above is

```{code-cell}
A = G.weighted_adjacency_matrix()
print(A)
```

If we are interested in the equilibrium probability of each state of the receptor model, it is sufficient to consider the rooted spanning tree

```{code-cell}
var('a b c kappa_b kappa_c')
T = DiGraph([[a,b,c],[(b,a),(c,b)]])
T.set_edge_label(b,a,kappa_b)
T.set_edge_label(c,b,kappa_c)
T.plot(figsize=4,edge_labels=True)
```

[To be completed]


# Scraps


```{code-cell}
B = T.weighted_adjacency_matrix()
print(B)
print(B**2)
print(B**3)
```

```{code-cell}
print(B)
print(B**2)
print(B**3)
```



```{code-cell}
G=graphs.PathGraph(3)
G.show(figsize=4)
```

```{code-cell}
G.relabel(dict({0: 'a', 1: 'b', 2: 'c'}))
G.show(figsize=4)
```

```{code-cell}
kappab = var("kappab", latex_name=r"\kappa_b")
kappac = var("kappac", latex_name=r"\kappa_c")
x = var("x", latex_name=r"x")
G.set_edge_label('a','b',kappab*x)
G.set_edge_label('b','c',kappab*kappac*x^2)
G.show(edge_labels=True,figsize=8)
f=kappab*kappac*x^2
f.show()
show(f)
print(f)
```

```{code-cell}
G=graphs.PathGraph(3).to_directed()
G.relabel({0:'R',1:'LR',2:'LLR'})
G.set_edge_label('R','LR','kap*L')
G.set_edge_label('LR','R','kam')
G.set_edge_label('LR','LLR','kbp*L')
G.set_edge_label('LLR','LR','kbm')
G.show(edge_labels=True,figsize=4,talk=True)
G.plot(edge_labels=True,figsize=4,talk=True)
```

## References 

```{bibliography}
:filter: docname in docnames
```


