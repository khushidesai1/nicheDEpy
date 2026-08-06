# nicheDEpy

`nicheDEpy` is a Python wrapper around the original
[NicheDE R package](https://github.com/kaishumason/NicheDE). It is designed for
Python workflows where spatial data are stored as `AnnData`, while preserving the
original NicheDE implementation for the statistical method.

The package does **not** reimplement NicheDE in Python. Instead, it:

1. extracts counts, coordinates, labels, and deconvolution weights from AnnData;
2. writes an R-readable working directory;
3. calls the upstream NicheDE R functions through `Rscript`;
4. returns standardized output paths and CSV tables to Python.

## Installation

Install the Python package:

```bash
pip install -e .
```

Install the R dependencies in an R session:

```r
install.packages(c("jsonlite", "Matrix"))
install.packages("devtools")
devtools::install_github("kaishumason/NicheDE")
```

If the upstream package is installed under a different GitHub owner in your
environment, install that version instead. `nicheDEpy` only requires that
`library(nicheDE)` works from `Rscript`.

## Quick Start

```python
import scanpy as sc
from nichedepy import run_pipeline

adata = sc.read_h5ad("spatial_data.h5ad")

result = run_pipeline(
    adata,
    labels_key="cell_type",
    spatial_key="spatial",
    sigma=[20, 200, 500],
    output_dir="outputs/nichede",
    num_cores=4,
)

print(result.object_path)
print(result.niche_de_object_path)
```

## Data Assumptions

`nicheDEpy` assumes:

- cells or spots are rows in `adata`;
- genes are columns in `adata`;
- spatial coordinates are in `adata.obsm[spatial_key]`;
- cell-type labels are in `adata.obs[labels_key]`;
- raw counts are either in `adata.layers[layer]` or `adata.X`;
- for single-cell data, deconvolution weights are generated as one-hot cell-type
  proportions from the labels.

For spot-level data, pass a deconvolution matrix with
`write_nichede_inputs(..., deconv=...)` before calling the lower-level wrappers.

## Main API

```python
from nichedepy import (
    write_nichede_inputs,
    create_library_matrix,
    create_niche_de_object,
    calculate_effective_niche,
    run_niche_de,
    get_niche_de_genes,
    run_pipeline,
)
```

The low-level functions mirror the user-facing R functions:

- `create_library_matrix()` wraps `CreateLibraryMatrix`.
- `create_niche_de_object()` wraps `CreateNicheDEObject`.
- `calculate_effective_niche()` wraps `CalculateEffectiveNiche`.
- `calculate_effective_niche(..., large_scale=True)` wraps
  `CalculateEffectiveNicheLargeScale`.
- `run_niche_de()` wraps `niche_DE`.
- `get_niche_de_genes()` wraps `get_niche_DE_genes`.
- `call_nichede_function()` provides an escape hatch for other exported NicheDE
  functions while keeping the original R implementation.

Additional named wrappers are also available for upstream helpers exported by
NicheDE, including `niche_de_markers`, `niche_de_no_parallel`, `filter_nde`,
`merge_objects`, `niche_lr_cell`, `niche_lr_spot`, `get_niche_de_pval_raw`,
`get_niche_de_pval_fisher`, `celltype_level`, `celltype_level_fisher`,
`gene_level`, `gene_level_fisher`, `check_colloc`, `contrast_post`, `nb_lik`,
`t_to_p`, and `ultosymmetric`.

## Example

See `examples/run_from_anndata.py` for an end-to-end template.

For users translating the official NicheDE R tutorial into Python, see
`docs/tutorial_tool_mapping.md`.
