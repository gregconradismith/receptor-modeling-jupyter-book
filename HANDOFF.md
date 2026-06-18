# Handoff

Current date: 2026-06-18

## Repository state

- Repository: `gregconradismith/receptor-modeling-jupyter-book`
- Local path: `/Users/gregconradismith/Library/CloudStorage/Dropbox/Main/Git/receptor-modeling-jupyter-book`
- Current branch: `upgrade`
- Remote tracking branch: `origin/upgrade`
- Latest pushed cleanup commit before this handoff: `47c00a2 Clean up book prose and archive drafts`

## What changed on `upgrade`

- Cleaned small prose issues and typos across active Jupyter Book pages.
- Standardized `SageMath` capitalization in touched prose.
- Renamed `hill_diagramatic_method` to `hill_diagrammatic_method` without aliases.
- Archived non-TOC draft/example material under `scraps/` instead of deleting it.
- Added `scraps/**` to `_config.yml` `exclude_patterns`.
- Left `README.md` as the only root-level source file not listed in `_toc.yml`.

## Important workflow note

The `deploy-book` workflow runs on:

- pushes to `main`
- pull requests targeting `main`
- manual `workflow_dispatch`

A push to `upgrade` alone does not trigger the workflow unless there is an open pull request from `upgrade` to `main`.

The workflow build job is pinned to `ubuntu-22.04` on this branch because `ubuntu-latest` currently resolves to Ubuntu 24.04, where `apt-get install sagemath` fails with no installation candidate.

## Local validation status

- Static TOC/source inventory was checked after archiving.
- Root-level non-TOC source inventory is now only `README`.
- Local `jupyter-book build .` was not run because `jupyter-book` is not installed in this shell.
- Local Sage also reported that the macOS SageMath app needs to be started once from the GUI before the `sage` command works.
- GitHub Actions is the best build validation path for this branch.

## Recommended next steps on the other computer

1. Fetch and check out the branch:

   ```bash
   git fetch origin
   git switch upgrade
   git pull
   ```

2. Open a pull request from `upgrade` to `main`, or manually dispatch `deploy-book` if you only want a build check.

3. Watch the GitHub Action:

   ```bash
   gh run list --workflow deploy.yml --limit 5
   gh run watch <run-id> --exit-status
   ```

4. If the PR build passes, merge to `main`; the deploy job should then publish to GitHub Pages.

Published site:

<https://gregconradismith.github.io/receptor-modeling-jupyter-book/>
