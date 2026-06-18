# Agent Instructions

This repository is a Jupyter Book about receptor modeling, allosteric
parameters, Markov chains, nonequilibrium steady states, cycle fluxes, and Wang
algebras.

Main files:

- `_config.yml` defines Jupyter Book settings, execution behavior, citations,
  repository links, and MyST extensions.
- `_toc.yml` defines the book structure.
- `README.md` documents local build and GitHub Pages publication.
- `.github/workflows/deploy.yml` is the canonical GitHub Pages build/deploy
  workflow for pushes to `main`.
- `requirements.txt` installs ordinary Python/Jupyter Book dependencies.
- `requirements-sage.txt` installs packages into SageMath's Python environment.
- Notebooks use the stable Jupyter kernel name `sagemath`.

Preserve the stable `sagemath` kernel name unless Greg explicitly changes the
workflow. Many notebooks are Sage-backed, so do not assume a plain Python kernel
is sufficient.

Be careful with notebooks:

- Avoid unnecessary execution-output churn.
- Keep notebook changes focused and review diffs before committing.
- Build or execute only the smallest relevant scope when practical; a full
  Jupyter Book build can be slow because it executes notebooks.

Useful checks:

```bash
git status --short --branch
git diff --check
jupyter-book build .
```

The GitHub Actions workflow is the canonical publication path. Local builds are
useful for diagnosing rendering and execution issues before pushing.

Do not commit local build artifacts or noise such as `_build/`, `.venv/`,
`.ipynb_checkpoints/`, `__pycache__/`, `.DS_Store`, editor temporaries, or
generated LaTeX artifacts.
