# coil-qc-code

This repository contains a reproducible Jupyter notebook for the manuscript titled "Multi-center benchmarking of cervical spinal cord RF coils at 7 T: A traveling spines study"

**Instructions to run the notebook locally:**

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

Download datasets
~~~
openneuro-py download --dataset ds005025 --target-dir data-human/
openneuro-py download --dataset ds005090 --target-dir data-phantom/
~~~


Run notebooks
~~~
jupyter notebook data_processing-human.ipynb
jupyter notebook data_processing-phantom.ipynb
~~~
