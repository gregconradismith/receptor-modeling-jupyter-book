# Codex Handoff

Date: 2026-06-20

Repo: `receptor-modeling-jupyter-book`

Branch: `main`

Current Git status at handoff update:

```bash
## main...origin/main
 D scraps/dated-notebooks/2023-08-20-example-lmfit.ipynb
 D scraps/dated-notebooks/2023-08-21-practicing-python.ipynb
 D scraps/dated-notebooks/2023-08-24-example-minimize-linear-regression.ipynb
 D scraps/dated-notebooks/2023-08-24-example-minimize-receptor-model-pade.ipynb
 D scraps/dated-notebooks/2023-08-25-plotting-hypercubes.ipynb
 D scraps/dated-notebooks/README.md
 D scraps/jupyter-book-examples/README.md
 D scraps/jupyter-book-examples/markdown-notebooks.md
 D scraps/jupyter-book-examples/markdown.md
 D scraps/legacy-mytools/EnumerateAllostericParameters.ipynb
 D scraps/legacy-mytools/MyStruggleWithDocstrings.ipynb
 D scraps/legacy-mytools/MyTools.ipynb
 D scraps/legacy-mytools/MyToolsExamplesEnumerateAllostericParameters.ipynb
 D scraps/legacy-mytools/MyToolsExamplesVertexAndEdgeLabels.ipynb
 D scraps/legacy-mytools/PlotBindingCurves.ipynb
 D scraps/legacy-mytools/README.md
 D scraps/legacy-mytools/enumerate_allosteric_parameters-2.ipynb
 D scraps/legacy-overviews/README.md
 D scraps/legacy-overviews/allosteric_parameters_overview.md
 D scraps/legacy-overviews/overview.md
 D scraps/legacy-overviews/wang_algebra.md
 D scraps/legacy-receptor-modeling/README.md
 D scraps/legacy-receptor-modeling/receptor_modeling.md
 D scraps/legacy-receptor-modeling/receptor_modeling_equilibrium_formalism.md
 D scraps/legacy-receptor-modeling/receptor_modeling_ligands.md
 D scraps/legacy-receptor-modeling/receptor_modeling_nonequilibrium.md
 D scraps/legacy-receptor-modeling/receptor_modeling_states_and_transitions.md
 D scraps/legacy-tool-examples/README.md
 D scraps/legacy-tool-examples/ReducedFromCartesianPower.ipynb
 D scraps/legacy-tool-examples/cycle_flux_derived_chain.ipynb
 D scraps/legacy-tool-examples/my_tools_examples_enumerate_allosteric_parameters.ipynb
 D scraps/legacy-tool-examples/my_tools_examples_vertex_and_edge_labels.ipynb
 D scraps/legacy-tool-examples/plot_binding_curves.ipynb
 D scraps/receptor_modeling_scraps.md
 M CODEX-HANDOFF.md
!! .DS_Store
```

The `scraps/**` deletions were already present when this handoff was updated.
Do not restore or commit them unless Greg confirms that is desired.

## Repository Role

This repository builds and publishes a Jupyter Book about receptor modeling,
including receptor-state graphs, ligands, thermodynamic constraints,
allosteric-parameter fitting, nonequilibrium steady states, cycle fluxes, and
Wang algebras.

## High-Value Context

- Read `AGENTS.md` before editing.
- `README.md` describes the local SageMath/Jupyter Book build setup.
- `_config.yml` uses cached notebook execution, `allow_errors: false`, and a
  600-second timeout.
- `_toc.yml` controls book order and chapter inclusion.
- Many notebooks require SageMath and should use a stable `sagemath` kernel.
- `.github/workflows/deploy.yml` builds on pushes and pull requests to `main`;
  deploy happens only for pushes to `main`.
- The workflow installs Jupyter Book dependencies with Python 3.12, installs
  SageMath on Ubuntu 22.04, registers the stable `sagemath` kernel name, builds
  the book, and deploys `_build/html` to GitHub Pages.
- A separate `upgrade` branch exists and contains an older `HANDOFF.md`; this
  root handoff was added directly on `main` per Greg's instruction to push
  `AGENTS.md` / `CODEX-HANDOFF.md` changes to `main`.

## Useful Commands

Create and activate the ordinary Python environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Install Sage-side dependencies:

```bash
sage -pip install -r requirements-sage.txt
```

Register the stable SageMath kernel name:

```bash
jupyter kernelspec list
jupyter kernelspec install --user --name sagemath /path/to/SageMath-kernel
```

Build the book:

```bash
jupyter-book build .
```

Check Git state:

```bash
git status --short --branch
git diff --check
```

## Notes For The Next Codex

- Preserve the `sagemath` kernel convention.
- Do not commit `_build/`, `.venv/`, notebook checkpoints, Python caches, local
  OS/editor files, or generated LaTeX artifacts.
- Avoid broad notebook-output churn.
- Prefer editing Markdown overview pages for prose/navigation changes and
  notebooks only when executable content needs to change.
- After adding or updating `AGENTS.md` / `CODEX-HANDOFF.md`, commit the scoped
  change on `main` and push to `origin/main`.

## Review Notes From 2026-06-20

Greg asked Codex to review the repository. No source fixes were made during
that review. Main findings:

- `_config.yml` enables only the `amsmath` MyST extension, while several pages
  still use `$...$` inline math. A non-executing HTML probe showed raw `$G$`,
  `$K_4$`, etc. Either enable `dollarmath` or convert remaining dollar math to
  `{math}` roles.
- Citation rendering needs cleanup. `receptor_modeling_overview.md` cites
  `Kenakin2018`, but the generated HTML showed an empty citation because no
  bibliography entry was emitted for that page. `wang_algebra_overview.md`
  cites `HammackSmith17`, while `references.bib` defines `HammackSmith2017`.
  Page-local bibliographies also duplicate keys, producing duplicate-citation
  warnings and cross-page citation targets.
- Published placeholders remain: `thermodynamic_constraints_and_allosteric_parameters.md`
  says "see XXX", and `further_reading.md` says "Add some further reading here."
- `equilibrium.md` links to `receptor_dimers_overview`, but that name is both a
  document stem and an explicit label, causing MyST ambiguous-reference warnings.
- `wang_algebra_directed_graph.ipynb` and `wang_algebra_undirected_graph.ipynb`
  jump from `#` to `###`, producing MyST header-level warnings.

Review commands run:

```bash
git diff --check
jupyter kernelspec list
jupyter-book build .
jupyter-book config sphinx .
sphinx-build -b html . /private/tmp/receptor-review-html -D nb_execution_mode=off
```

`git diff --check` passed. `jupyter-book build .` failed locally because this
machine did not have a stable `sagemath` kernelspec alias, even though the
GitHub Actions workflow registers one. The non-executing Sphinx build succeeded
and produced 18 warnings. The generated `conf.py` from `jupyter-book config
sphinx .` was deleted before the review ended.
