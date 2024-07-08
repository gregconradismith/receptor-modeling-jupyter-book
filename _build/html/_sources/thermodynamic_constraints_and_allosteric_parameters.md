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

(thermodynamic_constraints)=
# Thermodynamic constraints

The parameters for receptor models must statisfy thermodynamic contraints (see XXX for more discussion).  For equilibrium models, it suffices to consider a spanning tree of the state-transition model.  The equilibrium constant for each edge of this spanning tree can be freely chosen.  However, the equilibrium constant for the chords of the spanning tree are constrained by the fact that the free energy differences around any cycle must sum to zero.  That is, the product of the association constants must be unity.

For example, continue the ternery complex model of a G protein coupled receptor (GPCR).  The seven transmembreane rceptor (7TM) can bind to ligand as well as the alpha subunit of heterotrimeric GTP-binding protein.  These two binding processes result in a receptor model with four states: 


(ligands:ternary_complex_model)=

```{code-cell}
:label: my_ternary_complex_model_code
var('R LR RG LRG kap kam kbp kbm kcp kcm kdp kdm L G')
Gtcm = DiGraph({R: {LR:kap*L, RG:kcp*G}, LR: {R:kam, LRG:kbp*G}, LRG: {LR:kbm, RG:kdm}, RG: {R:kcm, LRG:kdp}})
pos = {R: (0, 0), RG: (2, 0), LR: (0, 2), LRG: (2, 2)}
Gtcm.plot(figsize=8,edge_labels=True,pos=pos,graph_border=True,vertex_size=1000)
```

The receptor model's state-transition diagram (above) has the topology of a symmetdfdic directed [path graph](example_graphs:path_graph) on 3 vertices. 

The transition `R` to `RL` occurs at rate `kap*L` where `L` is ligand concentration and `kap` is an association rate constant with physical dimensions of {math}`\mbox{time}^{-1} \mbox{conc}^{-1}` where {math}`\mbox{conc}` is concentration (i.e., number density, {math}`\mbox{amount}/\mbox{length}^3`). The transition `RL` to `L` that involes unbinding of ligand is unimolecular, i.e., physical dimensions of rate ({math}`\mbox{time}^{-1}`).  The interpretation of the edge weights `kbp*L` and `kbm` are similar. 

```{note}
The characters `kap` and `kam` stand for {math}`k_a^+` and {math}`k_a^-` (`p` for plus and `m` for minus). Similarly for `kbp` and `kbm`. The product of a bimolecular rate constant and ligand concentration `kap*L` stands for {math}`k_a^+ [{\rm L}]` where {math}`{\rm L}` is a chemical species and the brackets indicate the concentration (number density) of that species. 
```


