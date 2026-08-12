# Machine Learning and Topological Data Analysis for Clustering and Classification

This repository contains the replication codebase, Quarto analysis notebooks, and generated thesis outputs for the MSc Data & Computational Science Research Project at University College Dublin (UCD).

This project extends the Functional Topological Data Analysis (funTDA) framework by incorporating unsupervised clustering (K-Means, Hierarchical) and supervised machine learning classification (SVM, Random Forest, KNN) on persistence landscape topological features.

## Repository Structure

```text
funTDA/
|-- README.md                       # Main repository documentation
|-- LICENSE                         # MIT License (Catherine Higgins & Xiaoyu Qi)
|-- analysis/
|   |-- funTDA_ML.R                 # Primary R execution script for thesis replication
|   |-- notebooks_qmd/
|   |   |-- funTDA and ML.qmd       # Quarto document for automated report & output generation
|   |   `-- funTDA-and-ML.pdf       # Compiled PDF report
|   `-- supervisor/                 # Baseline funTDA framework scripts (reference only)
|       |-- funTDA.py               # Persistent homology computation (Python)
|       `-- funTDA.R                # Core FDA and persistence landscape script (R)
|-- data/
|   `-- GRN_adjacency_matrices/     # Raw H3N2 influenza gene regulatory network adjacency matrices
`-- output/
    |-- PH/                         # Extracted birth-death pair CSV datasets
    |-- R_figures/                  # Exported PDF diagnostic plots and panels
    `-- R_tables/                   # Exported LaTeX (.tex) summary and result tables
