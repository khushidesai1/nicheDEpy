from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from nichedepy import write_nichede_inputs


def test_write_nichede_inputs(tmp_path: Path) -> None:
    adata = ad.AnnData(
        X=np.array([[1, 0, 3], [0, 2, 1]], dtype=float),
        obs=pd.DataFrame({"cell_type": ["A", "B"]}, index=["c1", "c2"]),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
        obsm={"spatial": np.array([[0.0, 1.0], [2.0, 3.0]])},
    )

    output_dir = write_nichede_inputs(
        adata,
        tmp_path / "inputs",
        labels_key="cell_type",
    )

    assert (output_dir / "counts.mtx").exists()
    assert (output_dir / "coords.csv").exists()
    assert (output_dir / "labels.csv").exists()
    assert (output_dir / "deconv.csv").exists()

