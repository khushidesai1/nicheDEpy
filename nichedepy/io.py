"""AnnData input serialization for NicheDE."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


def _get_matrix(adata: ad.AnnData, layer: str | None) -> sparse.csr_matrix:
    """Return an AnnData expression matrix as CSR sparse matrix."""
    matrix: Any = adata.layers[layer] if layer is not None else adata.X
    if sparse.issparse(matrix):
        return matrix.tocsr()
    return sparse.csr_matrix(np.asarray(matrix))


def _get_spatial(adata: ad.AnnData, spatial_key: str) -> np.ndarray:
    """Return a two-column spatial coordinate array."""
    if spatial_key not in adata.obsm:
        raise KeyError(f"adata.obsm['{spatial_key}'] was not found.")
    coords = np.asarray(adata.obsm[spatial_key])
    if coords.ndim != 2 or coords.shape[1] < 2:
        raise ValueError("Spatial coordinates must have shape (n_obs, >=2).")
    return coords[:, :2]


def _one_hot_deconvolution(labels: pd.Series) -> pd.DataFrame:
    """Create one-hot deconvolution weights from cell-type labels."""
    deconv = pd.get_dummies(labels.astype(str), dtype=float)
    deconv.index = labels.index.astype(str)
    return deconv


def write_nichede_inputs(
    adata: ad.AnnData,
    output_dir: str | Path,
    labels_key: str,
    spatial_key: str = "spatial",
    layer: str | None = None,
    deconv: pd.DataFrame | None = None,
) -> Path:
    """Write AnnData fields to an R-readable NicheDE input directory.

    Args:
        adata: Input AnnData object with cells/spots as rows and genes as columns.
        output_dir: Directory where serialized inputs will be written.
        labels_key: Column in ``adata.obs`` containing cell-type labels.
        spatial_key: Key in ``adata.obsm`` containing spatial coordinates.
        layer: Optional AnnData layer to use for counts. If ``None``, ``adata.X``
            is used.
        deconv: Optional spot-by-cell-type deconvolution matrix. If omitted, a
            one-hot matrix is created from ``labels_key``.

    Returns:
        Path to the created input directory.
    """
    if labels_key not in adata.obs:
        raise KeyError(f"adata.obs['{labels_key}'] was not found.")

    output_path = Path(output_dir).expanduser().resolve()
    output_path.mkdir(parents=True, exist_ok=True)

    counts = _get_matrix(adata, layer)
    coords = _get_spatial(adata, spatial_key)
    obs_names = pd.Index(adata.obs_names.astype(str), name="cell_id")
    var_names = pd.Index(adata.var_names.astype(str), name="gene")
    labels = adata.obs[labels_key].astype(str).copy()
    labels.index = obs_names

    mmwrite(output_path / "counts.mtx", counts)
    pd.Series(obs_names, name="cell_id").to_csv(
        output_path / "cells.tsv", sep="\t", index=False
    )
    pd.Series(var_names, name="gene").to_csv(
        output_path / "genes.tsv", sep="\t", index=False
    )
    pd.DataFrame(coords, index=obs_names, columns=["x", "y"]).to_csv(
        output_path / "coords.csv"
    )
    pd.DataFrame({"cell_id": obs_names, "cell_type": labels.values}).to_csv(
        output_path / "labels.csv", index=False
    )

    if deconv is None:
        deconv = _one_hot_deconvolution(labels)
    else:
        deconv = deconv.copy()
        deconv.index = deconv.index.astype(str)
    deconv.to_csv(output_path / "deconv.csv")

    return output_path

