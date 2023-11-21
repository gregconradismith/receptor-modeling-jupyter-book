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
\def\b{\mathsf{b}}
\def\c{\mathsf{c}}
\def\kappab{\kappa_{\b}}
\def\kappac{\kappa_{\c}}
\def\kappabstar{\hat{\kappa}_{\b}}
\def\kappacstar{\hat{\kappa}_{\c}}
```

# The equilibrium formalism for receptor models 


Here we repeat our analysis of sequential binding using a notation that will prove helpful.

To find the equilibrium probability of each state of a receptor model, it is sufficient to consider the weighted rooted spanning tree that specifies the model.

```{code-cell}
var('a b c kappa_b kappa_c')
T = DiGraph([[a,b,c],[(b,a),(c,b)]])
T.set_edge_label(b,a,kappa_b)
T.set_edge_label(c,b,kappa_c)
T.plot(figsize=8,pos={a:(0,0),b:(2,0),c:(4,0)},edge_labels=True)
```

In the diagram above, {math}`\kappab` and {math}`\kappac` are dimensionless equilibrium constants, {math}`\kappabstar` and {math}`\kappacstar` are association constants with physical dimension of inverse concentration, and {math}`x` is ligand concentration.  

The edges of the spanning tree are directed backwards, i.e., the forward reaction direction is opposite of the direction of the arrow.  For example, the reaction labelled {math}`\kappab` has {math}`a` as reactant and {math}`b` as product; consequently, increasing {math}`\kappab` decreases the equilibrium probability (relative fraction) of state {math}`a` and increases the probability of state {math}`b`.

The three states are labelled so that the reactant comes before the product in dictionary order ({math}`a` to {math}`b` to {math}`c`).  The subscript of the equilibrium constants {math}`\kappab` and {math}`\kappac` are chosen to match the label of the _products_.


## Equilibrium analysis 


For the equilibrium receptor model above, the probability of state {math}`i` is given by
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



## References

```{bibliography}
:filter: docname in docnames
```


