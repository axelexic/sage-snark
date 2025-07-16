# sage-snark
Implementation of different SNARKs in SageMath. The goal of this repository to implement differnt proof systems for conceptual clarity over programming tricks for performance gains.

## Installation

You must have [SageMath](https://www.sagemath.org/) installed and available on your path to make use of these libraries. Run [./setup_python.sh](./setup_python.sh) to create a SageMath virtual environment, then run the following commands to install required python dependencies.

```shell
$ soure ./venv/bin/activate
$ pip install --upgrade pip && pip install -r requirements.txt
```

## Proof Systems

The following Jupyter notebooks provide detailed information about each of the supported proof systems:

1. [Spartan](./spartan.ipynb)

