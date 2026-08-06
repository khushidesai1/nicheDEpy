"""Configuration objects for NicheDE runs."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class NicheDEConfig:
    """Configuration for running NicheDE from AnnData.

    Attributes:
        labels_key: Column in ``adata.obs`` containing cell-type labels.
        spatial_key: Key in ``adata.obsm`` containing spatial coordinates.
        sigma: Kernel bandwidths passed to ``CreateNicheDEObject``.
        layer: Optional AnnData layer containing counts. If ``None``, ``adata.X``
            is used.
        output_dir: Directory used for serialized inputs and NicheDE outputs.
        num_cores: Number of cores passed to ``niche_DE``.
        rscript: Rscript executable.
        large_scale: Whether to use ``CalculateEffectiveNicheLargeScale``.
        batch_size: Batch size for large-scale effective niche calculation.
        cutoff: Cutoff for large-scale effective niche calculation.
    """

    labels_key: str
    spatial_key: str = "spatial"
    sigma: list[float] = field(default_factory=lambda: [20, 200, 500])
    layer: str | None = None
    output_dir: str | Path = "outputs/nichede"
    num_cores: int = 1
    rscript: str = "Rscript"
    large_scale: bool = False
    batch_size: int = 1000
    cutoff: float = 0.05

    def resolved_output_dir(self) -> Path:
        """Return the output directory as a resolved ``Path``."""
        return Path(self.output_dir).expanduser().resolve()

