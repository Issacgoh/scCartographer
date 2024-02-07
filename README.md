<img src="/resources/scCart_logo.webp?raw=true" alt="logo" width="100" height="100">

# scCartographer: Low-dimensional probabillistic transductive format transfer learning and mapping across latent spaces

## Introduction

`scCartographer` is a computational tool-set designed for the rapid calibrated, probabillistic transfer and mapping of categorical information across latent spaces derived from ldVAE (part of the scvi-tools suite). By leveraging an explainable set of latent spaces, it enables the training of a linear, Bayesian-optimized Elastic Net classifier for probabilistic category assignment between query and reference sets. This approach aims to provide insights into the linear decision-making processes of the model, enhancing the interpretability of complex biological data.

## Features

- **ldVAE Integration**: Utilizes ldVAE from the scvi-tools suite for high-quality latent space representations.
- **Bayesian-optimized Classifier**: Employs a Bayesian-optimized Elastic Net model for accurate and probabilistic category mapping.
- **Model Impact Assessment**: Includes a module for in-depth analysis of the factors influencing model decisions.
- **Comprehensive Analysis Tools**: Features ORA and Bayes DE modules for exploring both linear and non-linear data projections.

## Getting Started

To get started with `scCartographer`, install the package and follow the usage instructions for setting up your analysis environment.

```bash
pip install git+git@github.com:Issacgoh/scCartographer
```

## Usage

```python
from scCartographer import DataProcessor, Model, Analyzer

# Load and preprocess data
data = DataProcessor.load_data('path/to/data')
ldvae_space = DataProcessor.process_with_ldvae(data)

# Train the model and perform category mapping
model = Model.train(ldvae_space, data)
assignments = Model.predict_categories(model, data)

# Analyze model impact and perform ORA
impact = Analyzer.assess_model_impact(model, data)
ora_results = Analyzer.run_ora(data)
```

## Bayesian Optimization

The package employs Bayesian optimization to fine-tune the Elastic Net model, maximizing the efficiency of categorical assignments.

## References

`scCartographer` has been instrumental in several groundbreaking studies:

- "Yolk Sac Cell Atlas Reveals Multiorgan Functions During Human Early Development" in Science, DOI: 10.1126/science.add7564
- "Mapping the Developing Human Immune System Across Organs" as a preprint, DOI: 10.1101/2022.01.17.476665
- "Blood and Immune Development in Human Fetal Bone Marrow and Down Syndrome" in Nature, DOI: 10.1038/s41586-021-03929-x
- "Developmental Cell Programs are Co-opted in Inflammatory Skin Disease" in Science, DOI: 10.1126/science.aba6500

## Example outs
- scCartographer provides a framework for training Bayesian optimised label transfer models (Elasticnet) on the joint-latent representation of anchor states with integrated data for efficient mapping and harmonization of atlas data. scCartographer also allows user to study relative feature impact on modelling decisions and output probability 

![Alt text](/resources/cross_atlas_probabillistic_projection.png?raw=true "cross_atlas mapping and probabillistic projection")

- Training can be performed relative to any subset of the latent representation
![Alt text](/resources/general_model_self_projection.png?raw=true "single tissue mapping and probabillistic projection")

- scCartographer additionally provides modules to study relative feature impacts and associated gene expression pathways/modules on specific model decisions, allowing users to assess model relevance
![Alt text](/resources/Kera_positive_features.png?raw=true "Positive features")
![Alt text](/resources/Kera_negative_features.png?raw=true "Negative features")


- scCartographer provides modules to cluster and study predicted expression programs between shared/differential states based on weighted semantic similarity and shared genes within mixtures of any defined database (e.g GOBP, KEGG)
![Alt text](/resources/cluster_expression_prorgams_1.png?raw=true "Markov clusters semantic programs 1")
![Alt text](/resources/cluster_expression_prorgams_2.png?raw=true "Markov clusters semantic programs 1")

## Project Team

Led by Issac Goh at Newcastle University and the Sanger Institute. For more information about the team and their work, visit [Issac Goh's profile](https://haniffalab.com/team/issac-goh.html).

## Contact

For inquiries, please reach out to Issac Goh at ig7@sanger.ac.uk.

## Built With

- Scanpy
- scvi-tools
- ElasticNet from Scikit-learn
- Bayesian Optimization Libraries
