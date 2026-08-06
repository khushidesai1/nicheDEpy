# NicheDE Tutorial Tool Mapping

This guide maps the main tools used in the official
[NicheDE tutorial](https://kaishumason.github.io/NicheDE/articles/Tutorial.html)
to the corresponding `nicheDEpy` wrapper functions.

`nicheDEpy` keeps the original NicheDE R implementation intact. The Python
functions prepare AnnData inputs, call the R functions through `Rscript`, and
return saved outputs that can be loaded back into Python.

## Typical AnnData Workflow

For most single-cell spatial datasets stored as AnnData, use:

```python
from nichedepy import run_pipeline

result = run_pipeline(
    adata,
    labels_key="cell_type",
    spatial_key="spatial",
    sigma=[20, 200, 500],
    output_dir="outputs/nichede",
    num_cores=4,
)
```

This runs the same conceptual workflow as the tutorial:

1. build average cell-type expression profiles;
2. create a NicheDE object;
3. calculate effective niches;
4. run NicheDE differential expression.

## Tutorial Function Mapping

| Tutorial tool | `nicheDEpy` wrapper | What it connects to | What it does |
| --- | --- | --- | --- |
| `CreateLibraryMatrix()` | `create_library_matrix()` | AnnData expression matrix and cell-type labels serialized by `write_nichede_inputs()` | Computes the average expression profile for each cell type. This becomes the cell-type reference matrix used when creating the NicheDE object. |
| `CreateNicheDEObject()` | `create_niche_de_object()` | Counts, spatial coordinates, deconvolution/cell-type weights, sigma values, and the library matrix | Constructs the core NicheDE object containing expression, spatial coordinates, cell-type information, and kernel bandwidths. |
| `CalculateEffectiveNiche()` | `calculate_effective_niche()` | NicheDE object created by `create_niche_de_object()` | Computes effective niche covariates for each observation using the spatial kernel and cell-type composition. |
| `CalculateEffectiveNicheLargeScale()` | `calculate_effective_niche_large_scale()` or `calculate_effective_niche(..., large_scale=True)` | Same NicheDE object as above | Memory-conscious version of effective niche calculation for larger datasets. |
| `niche_DE()` | `run_niche_de()` | NicheDE object after effective niche calculation | Runs the main NicheDE differential expression model and stores the fitted results in the object. |
| `niche_DE_no_parallel()` | `niche_de_no_parallel()` | NicheDE object after effective niche calculation | Runs the NicheDE model without parallel execution. Useful for debugging or environments where parallel R execution is difficult. |
| `get_niche_DE_genes()` | `get_niche_de_genes()` | NicheDE result object after `run_niche_de()` | Extracts significant NicheDE genes for a selected index cell type, niche cell type, and test level. |
| `get_niche_DE_pval_raw()` | `get_niche_de_pval_raw()` | NicheDE result object | Extracts raw NicheDE p-values from fitted results. Use this when reproducing or customizing downstream result tables. |
| `get_niche_DE_pval_fisher()` | `get_niche_de_pval_fisher()` | NicheDE result object | Extracts Fisher-combined p-values from fitted results when using the corresponding NicheDE aggregation. |
| `Filter_NDE()` | `filter_nde()` | NicheDE result object | Filters NicheDE results by significance or result criteria defined by the upstream R function. |
| `niche_DE_markers()` | `niche_de_markers()` | NicheDE object/results | Runs marker-oriented NicheDE helper analysis exposed by the upstream package. |
| `niche_LR_cell()` | `niche_lr_cell()` | NicheDE object/results and ligand-receptor inputs expected by upstream NicheDE | Performs ligand-receptor analysis for cell-level NicheDE results. |
| `niche_LR_spot()` | `niche_lr_spot()` | NicheDE object/results and ligand-receptor inputs expected by upstream NicheDE | Performs ligand-receptor analysis for spot-level NicheDE results. |
| `MergeObjects()` | `merge_objects()` | Multiple NicheDE objects/results | Merges upstream NicheDE objects, matching the R helper behavior. |
| `celltype_level()` | `celltype_level()` | NicheDE test statistics/results | Aggregates or summarizes results at the cell-type level. |
| `celltype_level_fisher()` | `celltype_level_fisher()` | NicheDE test statistics/results | Cell-type-level summary using Fisher-style aggregation. |
| `gene_level()` | `gene_level()` | NicheDE test statistics/results | Aggregates or summarizes results at the gene level. |
| `gene_level_fisher()` | `gene_level_fisher()` | NicheDE test statistics/results | Gene-level summary using Fisher-style aggregation. |
| `check_colloc()` | `check_colloc()` | Spatial/cell-type information in the NicheDE workflow | Checks colocalization patterns using the upstream R helper. |
| `contrast_post()` | `contrast_post()` | NicheDE fitted results | Runs the upstream post-hoc contrast helper. |
| `nb_lik()` | `nb_lik()` | Model inputs expected by upstream NicheDE | Calls the upstream negative-binomial likelihood helper. Mostly useful for advanced users reproducing internals. |
| `T_to_p()` | `t_to_p()` | Test statistics | Converts test statistics to p-values using the upstream helper. |
| `ultosymmetric()` | `ultosymmetric()` | Matrix-like inputs | Converts an upper/lower triangular representation into a symmetric matrix using the upstream helper. |

## Input Preparation

The R tutorial starts from R matrices. In Python, the equivalent first step is:

```python
from nichedepy import write_nichede_inputs

input_dir = write_nichede_inputs(
    adata,
    output_dir="outputs/nichede/inputs",
    labels_key="cell_type",
    spatial_key="spatial",
    layer="counts",
)
```

This writes:

- `counts.mtx`: cells/spots by genes expression matrix;
- `coords.csv`: spatial coordinates;
- `labels.csv`: cell-type labels;
- `deconv.csv`: one-hot cell-type weights for single-cell data;
- `cells.tsv` and `genes.tsv`: matrix names.

For spot-level data with deconvolution weights, pass a spot-by-cell-type
`pandas.DataFrame` to `deconv`.

## Recommended Entry Points

Use `run_pipeline()` when you want the standard tutorial workflow from AnnData.

Use the lower-level functions when you want to inspect or modify intermediate
objects:

```python
input_dir = write_nichede_inputs(adata, "outputs/nichede/inputs", "cell_type")
library_path = create_library_matrix(input_dir)
object_path = create_niche_de_object(input_dir, sigma=[20, 200, 500])
effective_path = calculate_effective_niche(object_path)
results_path = run_niche_de(effective_path, num_cores=4)
genes = get_niche_de_genes(
    results_path,
    test_level="I",
    index="ReceiverCellType",
    niche="SenderCellType",
)
```

Use `call_nichede_function()` only when the upstream R package exposes a helper
that is not yet represented by a named Python wrapper.

