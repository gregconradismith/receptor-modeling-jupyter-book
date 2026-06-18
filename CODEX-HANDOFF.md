# Codex Handoff

Date: 2026-06-18

Repo: `receptor-modeling-jupyter-book`

Branch: `main`

Current Git status at handoff creation:

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
