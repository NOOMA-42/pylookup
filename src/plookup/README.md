# plookup

A Python implementation of [plookup](https://eprint.iacr.org/2020/315.pdf) with kzg proof

## Getting started

### 1. Install poetry

To get started, you'll need to have a Python version >= 3.8 and [`poetry`](https://python-poetry.org) installed:

`curl -sSL https://install.python-poetry.org | python3 -`.

### 2. Install dependencies

Run command in the root of the repository:

`poetry install`

This will install all the dependencies in a virtualenv.

### 3. Run

Run command from the root of the repository:

`poetry run python ./src/plookup/test.py`

This will take you through the workflow of setup, proof generation, and verification
