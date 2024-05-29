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
\def\bpit{\boldsymbol{\pi}^{\,T}}
\def\bzero{\boldsymbol{0}}
\def\be{\boldsymbol{e}}
```

(nonequilibrium)=
# Non-equilibrium steady states

The equilibrium formulation assumes detailed balance. However, when the state-transition diagram of a receptor model includes cycles, the steady-state probability distribution may not satisfy detailed balance.

Let us consider again the following state-transition diagram. 

```{code-cell} ipython3
var('a12, a21, a13, a31, a23, a32, a34, a43')
d = {1: {2:a12, 3:a13}, 2: {1:a21, 3:a23}, 3: {2:a32, 1:a31, 4:a34}, 4: {3:a43}};
G = DiGraph(d,weighted=True)
vertex_positions = {1: (0, 0), 2: (1, 1.41), 3: (2, 0), 4: (4,0)}
G.plot(figsize=8,edge_labels=True,pos=vertex_positions,graph_border=True)
```

## Generator matrix 

The generator matrix {math}`Q` for the Markov chain associated to {math}`G` can be constructed from the weighted adjacency matrix {math}`A`, as follows
```{code-cell}
A = G.weighted_adjacency_matrix(sparse=False)
Q = A - diagonal_matrix(sum(A.T))
show(Q)
print('The rank of Q is',Q.rank())
```

The generator of a Markov chain receptor model with {math}`n` states has rank {math}`n-1` due to conservation of probability. The four-state receptor model under consideration has a 4 states.  The {math}`4 \times 4` generator matrix {math}`Q` has rank {math}`4-1=3`.

```{code-cell}
print('The rank of Q is',Q.rank())
```

## The probability distribution solves a linear system of ordinary differential equations

The relationship between the adjacency matrix {math}`A` and generator matrix {math}`Q` can be written as
\begin{equation}
Q = A - \text{diag}(A\be)
\end{equation}
where {math}`\be` is a commensurate column vector of 1s. The following code confirms that {math}`Q \be = \bzero`, i.e., each row of {math}`Q` sums to zero.

```{code-cell}
e = matrix([1,1,1,1]).T
show(Q*e)
```

The matrix 
{math}`Q` is referred to as the generator matrix for the Markov chain because the probability distribution {math}`\bpit` solves 
```{math}
:label: mc_ode
d\bpit/dt = \bpit Q \, . 
```

```{note}
In the above expressions, the probability distribution {math}`\bpit` is a row vector that multiplies {math}`Q` on the left, while {math}`\be` is a column vector that multiplies {math}`Q` on the right.  If one prefers to represent the steady-state probability distribution as a column vector, one can write {math}`\be^T Q^T = \bzero^T` and {math}`d\bpi/dt = Q^T \bpi`.  
```

Setting the left side of {eq}`mc_ode` to zero, it is evident that steady-state probability distribution {math}`\bpit` solves {math}`\bpit Q = \bzero` subject to {math}`\bpit \be = 1`.  This expression is equivalent to {math}`\sum_i \pi_i = 1` (conservation of probability).


## Symbolic calcualtion of steady-state probability distribution


Using `Sagemath` we can symbolically solve {math}`\bpit Q = \bzero` subject to {math}`\bpit \be = 1`.  

To begin, we will unpack the four linear equations of {math}`\bpit Q = \bzero`, which is compact notation for four linear equations.

```{code-cell}
var('p1 p2 p3 p4')
p = vector([p1, p2, p3, p4])
pQ = p*Q
eq =[]
for lhs in pQ:
    show(lhs == 0)
    eq.append(lhs == 0)
```

As discussed above, the generator matrix {math}`Q` is rank 3.  Thus, the fourth equation (`a34*p3 - a43*p4 == 0`) is superfluous.  We replace this equation by the condition `p1+p2+p3+p4 == 1` (the solution should be a normalized probability distribution).

```{code-cell}
eq[-1] = p1+p2+p3+p4 == 1
for q in eq:
    show(q)
```

This system of four linear equations can now be solved to obtain the steady-state probability distribution. 

```{code-cell}
z = solve(eq,list(p))
for i in range(4):
    f = z[0][i].rhs()
    print('p%s' % (i+1),'=',f.expand().factor())
```

This solution can be written more compactly as follows.

```{math}
:label: Eq:P1234
p_1 & = \frac{z_1}{z_1+z_2+z_3+z_4}\\
p_2 & = \frac{z_2}{z_1+z_2+z_3+z_4}\\
p_3 & = \frac{z_3}{z_1+z_2+z_3+z_4}\\
p_4 & = \frac{z_4}{z_1+z_2+z_3+z_4}
```
where
\begin{align*}
z_1 & = (a_{21} a_{31} + a_{23} a_{31} + a_{21} a_{32} ) a_{43}\\
z_2 & = (a_{12} a_{31} + a_{12} a_{32} + a_{13} a_{32} ) a_{43}\\
z_3 & = (a_{13} a_{21} + a_{12} a_{23} + a_{13} a_{23} ) a_{43}\\
z_4 & = (a_{13} a_{21} + a_{12} a_{23} + a_{13} a_{23} ) a_{34} \, .
\end{align*}

In general, this probability distribution is a _non-equilibrium steady state_.  To see this, check to see if the distribution satisfies detailed balance, i.e., {math}`a_{ij} \pi_i = a_{ji} \pi_j`. Using the note above, we see that detailed balance implies {math}`a_{ij} z_i = a_{ji} z_j`, but {math}`a_{01} z_0 \neq a_{10} z_1` in general.


## [Komolgorov's criterion](https://en.wikipedia.org/wiki/Kolmogorov%27s_criterion) and equilibrium

If the product of rate constants around the cycle is the same in both directions, i.e., `a12*a23*a31=a13*a32*a21` [Komolgorov's criterion](https://en.wikipedia.org/wiki/Kolmogorov%27s_criterion) is satisfied.  In this case, the steady-state probability distribution is guaranteed to satisfy detailed balance.  The code below usings this condition to repace `a12` by `a13*a32*a21/(a23*a31)`.

```{code-cell} ipython3
:tags: ["hide-cell"]
def mysolve(p,Q):
    pQ = p*Q
    eq =[]
    for lhs in pQ:
       eq.append(lhs == 0)
    eq[-1] = p1+p2+p3+p4 == 1
    z = solve(eq,list(p))
    for i in range(4):
        f = z[0][i].rhs()
        print('p%s' % (i+1),'=',f.expand().factor())
```

The equilibrium steady-state probability distribution assuming the Komolgorov condition is
```{code-cell}
Q = Q.subs(a31=a13*a32*a21/(a12*a23))
mysolve(p,Q)
```
which can be written more compactly using {eq}`Eq:P1234` and
\begin{align*}
z_1 & = a_{21} a_{32} a_{43} &
z_2 & = a_{12} a_{32} a_{43} &
z_3 & = a_{12} a_{23} a_{43} &
z_4 & = a_{12} a_{23} a_{34} \, .
\end{align*}

This distribution satisfies detailed balance, i.e., {math}`a_{ij} \pi_i = a_{ji} \pi_j`.  As discussed previously, detailed balance implies that the steady-state probability distribution can be written in terms of equilibrium constants (as opposed to rate constants).  Define {math}`\kappa_{j}=a_{ij}/a_{ji}` for {math}`i < j` whenever vertex {math}`i` and {math}`j` are adjacent.
Dividing the numerator and denominator of each {math}`\pi_i` by {math}`a_{12}a_{23}a_{34}` yields the relative probabilities {math}`z_1 = 1`,   {math}`z_2 = \kappa_{2}`,  {math}`z_3 = \kappa_{2}\kappa_{3}`, and {math}`z_4 =  \kappa_{2} \kappa_{3} \kappa_{4}`.  Substituting these values into {eq}`Eq:P1234` gives the solution in terms of the equilibrium constants {math}`\kappa_{2}`, {math}`\kappa_{3}`, and {math}`\kappa_{4}`.

In the notation of the [equilibrium formulation](equilibrium), this probability distribution, 
\begin{equation}
 [ \pi_1  :  \pi_2 :  \pi_3 :  \pi_4] = [1 :\kappa_{2} : \kappa_{2} \kappa_{3} :  \kappa_{2} \kappa_{3} \kappa_{4} ]  \, ,
\end{equation}
corresponds to the following spanning tree rooted in state 1.
```{code-cell}
var('kappa_2, kappa_3, kappa_4')
d = {2: {1:kappa_2}, 3: {2:kappa_3}, 4: {3:kappa_4}};
G = DiGraph(d,weighted=True)
vertex_positions = {1: (0, 0), 2: (1, 1.41), 3: (2, 0), 4: (4,0)}
G.plot(figsize=8,edge_labels=True,pos=vertex_positions,graph_border=True)
```


