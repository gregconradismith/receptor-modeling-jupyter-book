# Receptor Modeling Jupyter Book

This repository contains a Jupyter Book about receptor modeling.

## Build locally

Install SageMath first. Most notebooks use the stable Jupyter kernel name
`sagemath`, so register your local SageMath kernel under that name before doing
an execution build.

There are two Python environments involved:

- the ordinary Python environment that runs `jupyter-book`
- the SageMath Python environment that executes Sage-backed notebooks

1. Create and activate the ordinary Python environment:

   ```bash
   python3 -m venv .venv
   source .venv/bin/activate
   python -m pip install -r requirements.txt
   ```

2. Install packages used by the SageMath notebooks into Sage's Python
   environment:

   ```bash
   sage -pip install -r requirements-sage.txt
   ```

   On macOS, if `sage` reports that the SageMath app cannot be found, start the
   SageMath app once from Finder and then retry the shell command.

3. Confirm that Jupyter can see a SageMath kernel:

   ```bash
   jupyter kernelspec list
   ```

   Look for a SageMath kernelspec path, for example:

   ```text
   sagemath-10.8    /usr/local/share/jupyter/kernels/SageMath-10.8
   ```

4. Register that SageMath kernelspec under the stable name `sagemath`:

   ```bash
   jupyter kernelspec install --user --name sagemath /usr/local/share/jupyter/kernels/SageMath-10.8
   ```

   Replace `/usr/local/share/jupyter/kernels/SageMath-10.8` with the SageMath
   path shown by `jupyter kernelspec list` on your Mac.

5. Confirm the stable kernel name is available:

   ```bash
   jupyter kernelspec list
   ```

   You should see an entry named `sagemath`.

6. Build the book:

   ```bash
   jupyter-book build .
   ```

The generated HTML site is written to `_build/html/index.html`.

## Publish with GitHub Pages

The repository includes a GitHub Actions workflow at
`.github/workflows/deploy.yml`. On pushes to `main`, it installs the book
dependencies, installs SageMath, registers the stable `sagemath` kernel name,
builds the book, and publishes `_build/html` to GitHub Pages.

The workflow is useful as the canonical publication path. Local builds are still
useful for rapid editing and for diagnosing notebook execution issues before
pushing.
