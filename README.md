# Multi-center benchmarking of cervical spinal cord RF coils at 7 T: A traveling spines study

[![example workflow](https://github.com/spinal-cord-7t/coil-qc-code/actions/workflows/run_notebooks.yml/badge.svg)](https://github.com/spinal-cord-7t/coil-qc-code/actions/workflows/run_notebooks.yml)
[![Dataset human](https://img.shields.io/badge/openneuro-human%20data-blue)](https://openneuro.org/datasets/ds005025)
[![Dataset phantom](https://img.shields.io/badge/openneuro-phantom%20data-yellow)](https://openneuro.org/datasets/ds005090)
[![Jupyter Book Badge](https://jupyterbook.org/badge.svg)](https://spinal-cord-7t.github.io/coil-qc-code)

This repository contains a reproducible Jupyter notebook for the manuscript titled "Multi-center benchmarking of cervical spinal cord RF coils at 7 T: A traveling spines study"

## See the notebook results

You can access the processed notebooks with outputs by clicking on the Jupyter Book badge 👆

## Run the notebook locally

### Pre-requesites

- (recommended) Linux or MacOS
- A virtual environment with Python 3.10.x installed
  - e.g. Using a [miniforge](https://github.com/conda-forge/miniforge) installation, `conda create -n coil-qc python=3.10` in a terminal then follow instructions for activation.
- Version 6.5 of the Spinal Cord Toolbox (SCT)
  - `git clone --depth 1 --single-branch --branch 6.5 https://github.com/spinalcordtoolbox/spinalcordtoolbox.git`
  - Follow the [SCT Installation instructions](https://spinalcordtoolbox.com/user_section/installation.html) for your operating system.

### Installation

Clone the [latest release](https://github.com/spinal-cord-7t/coil-qc-code/releases) of this repository, eg. for release `r20250319`:

~~~
git clone --depth 1 --single-branch --branch r20250319 https://github.com/spinal-cord-7t/coil-qc-code.git 
cd coil-qc-code
~~~

Install Python dependencies (ensure that the virtual environment is created, see pre-requesites section above)
~~~
pip install -r requirements.txt
~~~

Run the notebooks
~~~
jupyter notebook data_processing-human.ipynb
jupyter notebook data_processing-phantom.ipynb
~~~

Alternatively, you can execute the notebooks by building a Jupyter Book with `jupyter-book build .`:

~~~
jupyter-book build .
~~~

and then open the resulting `./_build/html/index.html` file in your web browser.