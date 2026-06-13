# HIPC Annotation Challenge: Infection Study Integration and Annotation as an example with scCRAFT+

This repository contains the code used to integrate infection study datasets and generate cell type annotations for the HIPC annotation challenge.

The main analysis is provided in:

```text
annotate_challenge.ipynb
```

The notebook performs four main steps:

1. Load infection study `.h5ad` files.
2. Standardize feature names across studies.
3. Integrate and annotate cells using `scCRAFTplus`.
4. Map predicted labels to the provided ontology cell type labels and export annotation CSV files.

---

## Repository Structure

```text
.
├── annotate_challenge.ipynb
└── README.md
```

---

## Input Data


Each study should be stored as:

```text
recipe_data/
├── infection_study_01/
│   └── infection_study_01_processed.h5ad
├── infection_study_02/
│   └── infection_study_02_processed.h5ad
├── infection_study_03/
│   └── infection_study_03_processed.h5ad
├── infection_study_04/
│   └── infection_study_04_processed.h5ad
└── infection_study_06/
    └── infection_study_06_processed.h5ad
```


## Required Python Packages

The notebook requires the following Python packages:

```text
numpy
pandas
scanpy
scipy
scikit-learn
torch
scCRAFTplus
```

The main imports used in the notebook are:

```python
import os
import copy
import random
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from scipy.stats import pearsonr, spearmanr
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader

import scCRAFTplus
from scCRAFTplus.annot import *
from scCRAFTplus.model import *
```

---

## How to Run

Open the notebook:

```bash
jupyter notebook annotate_challenge.ipynb
```

or:

```bash
jupyter lab annotate_challenge.ipynb
```

Then run the cells from top to bottom.

Before running, check and update the following paths if needed:

```python
base_dir = "~/immunespace.org/hipc-sc-comp-data/team01/recipe_data"
```

and:

```python
out_dir = "~/results"
```

---

## Method Overview

### 1. Load infection studies

The notebook loads the selected infection studies:

```python
study_ids = ["01", "02", "03", "04", "06"]
```

For each study, the corresponding `.h5ad` file is loaded from:

```python
infection_study_{sid}/infection_study_{sid}_processed.h5ad
```

Each cell barcode is prefixed with its study name to make cell IDs unique after merging:

```python
adata.obs_names = [f"{study_name}_{x}" for x in adata.obs_names.astype(str)]
```

A study label is also added:

```python
adata.obs["study"] = study_name
```

---

### 2. Standardize feature names

The notebook includes a helper function:

```python
standardize_var_names(adata, study_name)
```

This function checks whether feature names are valid and harmonizes feature naming across studies.

The function handles the following cases:

* If `var_names` are numeric values such as `0`, `1`, `2`, etc., the notebook raises an error because those are not valid gene names.
* If `var_names` are Ensembl IDs and a `feature_name` column is available in `.var`, the notebook replaces `var_names` with `adata.var["feature_name"]`.
* If `var_names` already look like gene symbols, they are kept unchanged.
* Empty or invalid feature names are removed.
* Duplicated gene names are made unique using:

```python
adata.var_names_make_unique()
```

This step is important because the studies are merged by shared features.

---

### 3. Merge studies by shared genes

After loading and standardizing all selected studies, the notebook merges them using:

```python
adata_merged_infection = sc.concat(
    adata_list,
    join="inner",
    axis=0,
    label="study_concat",
    keys=study_keys,
    index_unique=None
)
```

The option:

```python
join="inner"
```

keeps only genes shared across all selected studies.

The merged object is then assigned to:

```python
adata = adata_merged_infection
```

---

### 4. Preprocess data

The notebook performs basic preprocessing before annotation and integration:

```python
adata.raw = adata
sc.pp.filter_cells(adata, min_genes=300)
sc.pp.filter_genes(adata, min_cells=5)
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
```

Then marker ranking and marker-based processing are performed using the PBMC level 3 marker list:

```python
ranks = rank_genes(adata, markers_PBMC_l3)
process_cell_types_with_ranks(adata, markers_PBMC_l3, ranks)
```

Highly variable genes are selected with study-aware batching:

```python
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="study")
adata = adata[:, adata.var["highly_variable"]]
```

---

### 5. Train scCRAFTplus model

The notebook clusters the data and trains the integration model:

```python
multi_resolution_cluster(adata, resolution1=0.5, method="Louvain")

VAE, C = train_integration_model(
    adata,
    batch_key="sample_id",
    z_dim=50,
    epochs=150,
    warmup_epoch=50
)
```

The trained model is then used to generate embeddings and predicted labels:

```python
obtain_embeddings(
    adata,
    VAE,
    C,
    cell_types_markers=markers_PBMC_l3
)
```

UMAP visualization is generated with:

```python
sc.tl.umap(adata, min_dist=0.5)
sc.pl.umap(adata, color=["sample_id", "Predict_label"], frameon=False, ncols=1)
sc.pl.umap(adata, color=["study", "Predict_label"], frameon=False, ncols=1)
```

---

### 6. Save annotated AnnData object

Before saving, the notebook converts the `age` column to string to avoid HDF5 writing issues caused by mixed data types:

```python
adata.obs["age"] = adata.obs["age"].astype(str)
```

The annotated AnnData object is saved as:

```python
adata.write_h5ad(
    "~/infection_filter_scCRAFT.h5ad"
)
```

---

### 7. Map predicted labels to ontology labels

The notebook maps `adata.obs["Predict_label"]` to the ontology-compatible `cell_type` column using a manual dictionary:

```python
adata.obs["cell_type"] = adata.obs["Predict_label"].map(predict_to_celltype)
```

The mapping collapses detailed predicted labels into broader ontology labels.

Examples include:

```python
"B_naive_kappa" -> "Naive B Cell"
"B_memory_kappa" -> "Memory B Cell"
"CD14_Mono" -> "Classical Monocyte"
"CD16_Mono" -> "Non-Classical Monocyte"
"CD4_TCM_1" -> "CD4 Naive / T Central Memory"
"CD8_TEM_1" -> "CD8 Cytotoxic / T Effector Memory"
"NK_1" -> "NK Cell"
"cDC1" -> "Conventional DC 1"
"cDC2_1" -> "Conventional DC 2"
"pDC" -> "Plasmacytoid DC"
"Platelet" -> "Platelet"
"Eryth" -> "RBC"
"HSPC" -> "HSC"
```

The notebook checks for unmapped predicted labels:

```python
unmapped = sorted(
    adata.obs.loc[adata.obs["cell_type"].isna(), "Predict_label"]
    .dropna()
    .unique()
)

print("Number of unmapped Predict_label values:", len(unmapped))
print("Unmapped labels:")
print(unmapped)
```

---

## Output Files

The notebook writes one annotation CSV file per study to:

```text
~/results
```

Each output file is named:

```text
infection_study_XX_annotation.csv
```

For example:

```text
infection_study_01_annotation.csv
infection_study_02_annotation.csv
infection_study_03_annotation.csv
infection_study_04_annotation.csv
infection_study_06_annotation.csv
```

Each CSV file contains one column:

```text
cell_type
```

The cell barcodes are stored as the row index.

Example output format:

```text
,cell_type
AAACCCAAGGGCAATC-1,Naive B Cell
AAAGGATTCTTCGACC-1,Memory B Cell
```

The study prefix is removed from the cell barcode before saving. For example:

```text
infection_study_01_AAACCCAAGGGCAATC-1
```

is saved as:

```text
AAACCCAAGGGCAATC-1
```

---

## Notes

The notebook assumes that the `scCRAFTplus` package and the PBMC marker list `markers_PBMC_l3` are available in the Python environment.

The notebook also assumes that the input `.h5ad` files contain the expected metadata columns, including:

```text
sample_id
study
age
```

If a dataset has inconsistent feature names, the function `standardize_var_names()` should be checked before integration.

If an output path is different from the default project directory, update:

```python
out_dir = "~/results"
```

before running the final export cell.

---

## Main Result

The final submitted annotation files are generated from the `cell_type` column, which is produced by mapping the model-predicted `Predict_label` values to the provided ontology labels.

The final per-study annotation files are saved as:

```text
infection_study_XX_annotation.csv
```

with one annotation column:

```text
cell_type
```
