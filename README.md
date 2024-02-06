# scCartographer
### Multi-modal Probabillistic Ensemble Transductive Transfer Learning

This repository contains a prototype classifier developed for mapping large sc-multiomic datasets using pre-integrated low-dimensional latent embeddings as input. If given input is linear, a weighted score for feature importance of predicted classes may be derived. Else, we can estimate this using a bayesian DE approach.

- Additionally includes a GCN-type attention-based framework for learning underlying graph structures in sc MD data

## Installation

You can install this package by cloning the GitHub repository and running `python setup.py install`.

## Usage

The main parts of the package are the data processing, modeling, and post-processing modules. Here's a basic example of how to use the package:

```python
from ** import data_processing, modeling, post_processing, feature_estimator
