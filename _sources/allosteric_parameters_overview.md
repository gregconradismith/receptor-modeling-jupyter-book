
# Overview

Quantitative pharmacologists construct Markov chain models to give insight into the relationship between ligand concentration and the fraction of cell surface receptors in each of several molecular conformations. Pharmacologists use these stochastic models to understand the action of natural ligands and drugs on receptor-mediated cell responses. When receptors function as two or more similar protein subunits working in concert (i.e., homodimers or oligomers), receptor models must

* account for symmetry,
* satisfy thermodynamic constraints, and
* properly account for subunit interactions (allostery) mediated by conformational coupling.

The modeling framework that satisfies these three requirements will be explicated in the context of models of G protein-coupled receptors (GPCRs), such as metabotropic glutamate receptors, that function as multi-molecule signaling complexes.  For equilibrium models of receptor dimers, this approach facilitates the inference of a parsimonious subset of allosteric interactions leading to conformational coupling and dependence of receptor subunits.

This manual presents this methodology and introduces helpful notation that was developed as part of this work. *En passant*, we distinguish two ways that thermodynamic constraints and allosteric parameters arise in receptor models:

* when the state-transition graph of a receptor model includes cycles, and

* as an emergent property of conformational coupling of receptor oligomers.


## Prerequisites

This manual presumes an understanding of receptor modeling and mathematical concepts familiar to the mainstream pharmacological community; see this: {cite}`Kenakin2018` for review.

## References 

```{bibliography}
:filter: docname in docnames
```
