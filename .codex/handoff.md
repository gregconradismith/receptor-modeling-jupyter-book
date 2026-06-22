# Codex Handoff

Date: 2026-06-22

Repo: `receptor-modeling-jupyter-book`

Branch: `main`

<!-- codex-transfer-snapshot:start -->
## 2026-06-22 Computer Transfer Snapshot

- Checked on 2026-06-22 from `/Users/greg/Git` before moving computers.
- Ran `git fetch --all --prune`; `main` is tracking `origin/main` unless this status says otherwise.
- Origin: `git@github.com:gregconradismith/receptor-modeling-jupyter-book.git`
- Latest commit at refresh time: `b9f59f8 2026-06-21 15:45:32 -0400 Refresh Codex handoff for computer migration`
- On the next machine, read `AGENTS.md` first, then this handoff.
- The working tree was clean before this handoff refresh; after committing the refresh, `git status --short --branch` should again show only the branch line.

Status before this handoff edit:

```bash
## main...origin/main
```
<!-- codex-transfer-snapshot:end -->

Current Git status after the 2026-06-21 migration readiness fetch and before this handoff edit:

```bash
## main...origin/main
```

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
  `AGENTS.md` / `.codex/handoff.md` changes to `main`.

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

## Migration Readiness Snapshot

- Checked on 2026-06-21 before moving computers.
- Non-interactive `git fetch --all --prune` completed successfully.
- Root `README.md` points to `.codex/handoff.md` when a root README exists.

Pre-edit Git state after fetch:

```bash
## main...origin/main
```

## Notes For The Next Codex

- Preserve the `sagemath` kernel convention.
- Do not commit `_build/`, `.venv/`, notebook checkpoints, Python caches, local
  OS/editor files, or generated LaTeX artifacts.
- Avoid broad notebook-output churn.
- Prefer editing Markdown overview pages for prose/navigation changes and
  notebooks only when executable content needs to change.
- After adding or updating `AGENTS.md` / `.codex/handoff.md`, commit the scoped
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
