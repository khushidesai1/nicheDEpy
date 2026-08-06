"""Run NicheDE from an AnnData object.

Usage:
    python examples/run_from_anndata.py path/to/data.h5ad cell_type outputs/nichede
"""

from __future__ import annotations

import sys

import anndata as ad

from nichedepy import run_pipeline


def main() -> None:
    """Run the standard NicheDE pipeline."""
    adata_path = sys.argv[1]
    labels_key = sys.argv[2]
    output_dir = sys.argv[3]

    adata = ad.read_h5ad(adata_path)
    result = run_pipeline(
        adata,
        labels_key=labels_key,
        spatial_key="spatial",
        sigma=[20, 200, 500],
        output_dir=output_dir,
        num_cores=4,
        large_scale=True,
    )
    print(result)


if __name__ == "__main__":
    main()

