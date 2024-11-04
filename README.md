# Multi-center benchmarking of cervical spinal cord RF coils at 7 T: A traveling spines study

[![example workflow](https://github.com/spinal-cord-7t/coil-qc-code/actions/workflows/run_notebooks.yml/badge.svg)](https://github.com/spinal-cord-7t/coil-qc-code/actions/workflows/run_notebooks.yml)
[![Dataset human](https://img.shields.io/badge/openneuro-human%20data-blue)](https://openneuro.org/datasets/ds005025)
[![Dataset phantom](https://img.shields.io/badge/openneuro-phantom%20data-yellow)](https://openneuro.org/datasets/ds005090)
[![Jupyter Book Badge](https://jupyterbook.org/badge.svg)](https://spinal-cord-7t.github.io/coil-qc-code)

This repository contains a reproducible Jupyter notebook for the manuscript titled "Multi-center benchmarking of cervical spinal cord RF coils at 7 T: A traveling spines study"

## See the notebook results

You can access the processed notebooks with outputs by clicking on the Jupyter Book badge 👆

## Run the notebook locally

Install [Spinal Cord Toolbox](https://spinalcordtoolbox.com/user_section/installation.html)

Clone this repository
~~~
git clone https://github.com/spinal-cord-7t/coil-qc-code.git 
cd coil-qc-code
~~~

Install Python dependencies (assuming Python is already installed)
~~~
pip install -r requirements.txt
~~~

Run notebooks
~~~
jupyter notebook data_processing-human.ipynb
jupyter notebook data_processing-phantom.ipynb
~~~
