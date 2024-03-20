
(receptor_dimers_overview)=
# Overview

Quantitative pharmacologists construct Markov chain models to give insight into the relationship between ligand concentration and the fraction of cell surface receptors in each of several molecular conformations. Pharmacologists use these stochastic models to understand the action of natural ligands and drugs on receptor-mediated cell responses. 

When receptors function as two or more similar protein subunits working in concert (i.e., homodimers or oligomers), a composite receptor model must


* account for symmetry,
* satisfy thermodynamic constraints, and
* properly account for subunit interactions (allostery) mediated by conformational coupling.

This section illustrates an approach to modeling conformational coupling of receptor dimers that was introduced in {cite:p}`HammackSmith2017` and {cite:p}`ConradiSmith2020`.

* Gregory Douglas Conradi Smith, **Allostery in oligomeric receptor models**, *Mathematical Medicine and Biology: A Journal of the IMA*, 37(3):313-333, 2020. [doi: 10.1093/imammb/dqz016](https://doi.org/10.1093/imammb/dqz016)

* Richard H. Hammack and Gregory D. Smith, **Cycle bases of reduced powers of graphs**, *ARS Mathematica Contemporanea*, 12(1):183–203, 2017. [doi: 10.26493/1855-3974.856.4d2](https://doi.org/10.26493/1855-3974.856.4d2)


For equilibrium models of receptor dimers, this approach facilitates the inference of a parsimonious subset of allosteric interactions leading to conformational coupling and dependence of receptor subunits.

Our first step is to introduce the concept of a [Cartesian product graph](https://en.wikipedia.org/wiki/Cartesian_product_of_graphs), because this clarifies the structure of the state-transition diagram of a homodimer composed of two identical (but not necessarily independent) monomeric subunits.

## References 

```{bibliography}
:filter: docname in docnames
```
