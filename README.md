# Receptor Modeling Jupyter Book

This repository contains a Jupyter Book about receptor modeling.

## Build on macOS

Install SageMath first. The notebooks use the stable Jupyter kernel name
`sagemath`, so register your local SageMath kernel under that name before
building.

1. Create and activate a Python environment:

   ```bash
   python3 -m venv .venv
   source .venv/bin/activate
   python -m pip install -r requirements.txt
   ```

2. Confirm that Jupyter can see a SageMath kernel:

   ```bash
   jupyter kernelspec list
   ```

   Look for a SageMath kernelspec path, for example:

   ```text
   sagemath-10.8    /usr/local/share/jupyter/kernels/SageMath-10.8
   ```

3. Register that SageMath kernelspec under the stable name `sagemath`:

   ```bash
   jupyter kernelspec install --user --name sagemath /usr/local/share/jupyter/kernels/SageMath-10.8
   ```

   Replace `/usr/local/share/jupyter/kernels/SageMath-10.8` with the SageMath
   path shown by `jupyter kernelspec list` on your Mac.

4. Confirm the stable kernel name is available:

   ```bash
   jupyter kernelspec list
   ```

   You should see an entry named `sagemath`.

5. Build the book:

   ```bash
   jupyter-book build .
   ```

The generated HTML site is written to `_build/html/index.html`.
