"""Public Python API for NicheDE."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import anndata as ad
import pandas as pd

from nichedepy.config import NicheDEConfig
from nichedepy.io import write_nichede_inputs
from nichedepy.runner import run_r_bridge


@dataclass
class NicheDEResult:
    """Paths produced by a NicheDE pipeline run."""

    input_dir: Path
    library_matrix_path: Path
    object_path: Path
    effective_niche_object_path: Path
    niche_de_object_path: Path


def create_library_matrix(
    input_dir: str | Path,
    output_path: str | Path | None = None,
    *,
    rscript: str = "Rscript",
) -> Path:
    """Create an average expression profile matrix with NicheDE.

    This wraps the R function ``CreateLibraryMatrix``.

    Args:
        input_dir: Directory created by ``write_nichede_inputs``.
        output_path: Optional CSV path for the library matrix.
        rscript: Rscript executable.

    Returns:
        Path to the written library matrix CSV.
    """
    input_path = Path(input_dir).expanduser().resolve()
    output = Path(output_path or input_path / "library_matrix.csv").resolve()
    run_r_bridge(
        "create_library_matrix",
        {"input_dir": str(input_path), "output_path": str(output)},
        rscript=rscript,
    )
    return output


def create_niche_de_object(
    input_dir: str | Path,
    sigma: list[float],
    output_path: str | Path | None = None,
    *,
    library_matrix_path: str | Path | None = None,
    rscript: str = "Rscript",
) -> Path:
    """Create a NicheDE object from serialized AnnData inputs.

    This wraps the R function ``CreateNicheDEObject``.

    Args:
        input_dir: Directory created by ``write_nichede_inputs``.
        sigma: Kernel bandwidths passed to NicheDE.
        output_path: Optional RDS path for the object.
        library_matrix_path: Optional existing library matrix CSV.
        rscript: Rscript executable.

    Returns:
        Path to the written NicheDE object RDS.
    """
    input_path = Path(input_dir).expanduser().resolve()
    output = Path(output_path or input_path / "nichede_object.rds").resolve()
    library_path = (
        Path(library_matrix_path).expanduser().resolve()
        if library_matrix_path is not None
        else input_path / "library_matrix.csv"
    )
    run_r_bridge(
        "create_object",
        {
            "input_dir": str(input_path),
            "sigma": sigma,
            "library_matrix_path": str(library_path),
            "output_path": str(output),
        },
        rscript=rscript,
    )
    return output


def calculate_effective_niche(
    object_path: str | Path,
    output_path: str | Path | None = None,
    *,
    large_scale: bool = False,
    batch_size: int = 1000,
    cutoff: float = 0.05,
    rscript: str = "Rscript",
) -> Path:
    """Calculate effective niche values for a NicheDE object.

    This wraps ``CalculateEffectiveNiche`` or
    ``CalculateEffectiveNicheLargeScale``.

    Args:
        object_path: Input NicheDE object RDS.
        output_path: Optional output RDS path.
        large_scale: Use the large-scale implementation.
        batch_size: Batch size for the large-scale implementation.
        cutoff: Cutoff for the large-scale implementation.
        rscript: Rscript executable.

    Returns:
        Path to the updated NicheDE object RDS.
    """
    object_file = Path(object_path).expanduser().resolve()
    output = Path(output_path or object_file.with_name("nichede_effective.rds"))
    run_r_bridge(
        "calculate_effective_niche",
        {
            "object_path": str(object_file),
            "output_path": str(output),
            "large_scale": large_scale,
            "batch_size": batch_size,
            "cutoff": cutoff,
        },
        rscript=rscript,
    )
    return output


def calculate_effective_niche_large_scale(
    object_path: str | Path,
    output_path: str | Path | None = None,
    *,
    batch_size: int = 1000,
    cutoff: float = 0.05,
    rscript: str = "Rscript",
) -> Path:
    """Calculate effective niche values with the large-scale R implementation.

    This is a convenience wrapper around ``CalculateEffectiveNicheLargeScale``.
    """
    return calculate_effective_niche(
        object_path,
        output_path,
        large_scale=True,
        batch_size=batch_size,
        cutoff=cutoff,
        rscript=rscript,
    )


def run_niche_de(
    object_path: str | Path,
    output_path: str | Path | None = None,
    *,
    num_cores: int = 1,
    outfile: str = "",
    C: float = 150,
    M: float = 10,
    gamma: float = 0.8,
    print_progress: bool = True,
    Int: bool = True,
    batch: bool = True,
    self_EN: bool = False,
    G: float = 1,
    rscript: str = "Rscript",
) -> Path:
    """Run NicheDE differential expression.

    This wraps the R function ``niche_DE``.

    Args:
        object_path: Input NicheDE object RDS after effective niche calculation.
        output_path: Optional output RDS path.
        num_cores: Number of cores for NicheDE.
        outfile: NicheDE output file prefix passed to R.
        C: NicheDE ``C`` parameter.
        M: NicheDE ``M`` parameter.
        gamma: NicheDE ``gamma`` parameter.
        print_progress: Whether NicheDE prints progress.
        Int: Whether to calculate interaction-level statistics.
        batch: Whether to use batch mode in NicheDE.
        self_EN: Whether to include self effective niches.
        G: NicheDE ``G`` parameter.
        rscript: Rscript executable.

    Returns:
        Path to the NicheDE object RDS after differential expression.
    """
    object_file = Path(object_path).expanduser().resolve()
    output = Path(output_path or object_file.with_name("nichede_results.rds"))
    run_r_bridge(
        "niche_de",
        {
            "object_path": str(object_file),
            "output_path": str(output),
            "num_cores": num_cores,
            "outfile": outfile,
            "C": C,
            "M": M,
            "gamma": gamma,
            "print": print_progress,
            "Int": Int,
            "batch": batch,
            "self_EN": self_EN,
            "G": G,
        },
        rscript=rscript,
    )
    return output


def niche_de_no_parallel(
    object_path: str | Path,
    output_path: str | Path,
    *,
    args: list[Any] | None = None,
    kwargs: dict[str, Any] | None = None,
    rscript: str = "Rscript",
) -> Path | None:
    """Wrap the upstream ``niche_DE_no_parallel`` helper."""
    return call_nichede_function(
        "niche_DE_no_parallel", object_path, output_path, args=args, kwargs=kwargs, rscript=rscript
    )


def niche_de_markers(
    object_path: str | Path,
    output_path: str | Path,
    *,
    args: list[Any] | None = None,
    kwargs: dict[str, Any] | None = None,
    rscript: str = "Rscript",
) -> Path | None:
    """Wrap the upstream ``niche_DE_markers`` helper."""
    return call_nichede_function(
        "niche_DE_markers", object_path, output_path, args=args, kwargs=kwargs, rscript=rscript
    )


def get_niche_de_genes(
    object_path: str | Path,
    test_level: str,
    index: str,
    niche: str,
    *,
    positive: bool = True,
    alpha: float = 0.05,
    output_path: str | Path | None = None,
    rscript: str = "Rscript",
) -> pd.DataFrame:
    """Return NicheDE genes for an index and niche cell type.

    This wraps the R function ``get_niche_DE_genes``.

    Args:
        object_path: NicheDE result object RDS.
        test_level: NicheDE test level, e.g. ``"I"``.
        index: Index cell type.
        niche: Niche cell type.
        positive: Return positive or negative NicheDE genes.
        alpha: FDR threshold.
        output_path: Optional CSV output path.
        rscript: Rscript executable.

    Returns:
        DataFrame containing the returned NicheDE genes.
    """
    object_file = Path(object_path).expanduser().resolve()
    output = Path(
        output_path
        or object_file.with_name(f"nichede_genes_{index}_{niche}_{test_level}.csv")
    )
    run_r_bridge(
        "get_niche_de_genes",
        {
            "object_path": str(object_file),
            "output_path": str(output),
            "test_level": test_level,
            "index": index,
            "niche": niche,
            "positive": positive,
            "alpha": alpha,
        },
        rscript=rscript,
    )
    return pd.read_csv(output)


def get_niche_de_pval_raw(
    object_path: str | Path,
    output_path: str | Path,
    *,
    args: list[Any] | None = None,
    kwargs: dict[str, Any] | None = None,
    rscript: str = "Rscript",
) -> Path | None:
    """Wrap the upstream ``get_niche_DE_pval_raw`` helper."""
    return call_nichede_function(
        "get_niche_DE_pval_raw", object_path, output_path, args=args, kwargs=kwargs, rscript=rscript
    )


def get_niche_de_pval_fisher(
    object_path: str | Path,
    output_path: str | Path,
    *,
    args: list[Any] | None = None,
    kwargs: dict[str, Any] | None = None,
    rscript: str = "Rscript",
) -> Path | None:
    """Wrap the upstream ``get_niche_DE_pval_fisher`` helper."""
    return call_nichede_function(
        "get_niche_DE_pval_fisher",
        object_path,
        output_path,
        args=args,
        kwargs=kwargs,
        rscript=rscript,
    )


def filter_nde(
    object_path: str | Path,
    output_path: str | Path,
    *,
    args: list[Any] | None = None,
    kwargs: dict[str, Any] | None = None,
    rscript: str = "Rscript",
) -> Path | None:
    """Wrap the upstream ``Filter_NDE`` helper."""
    return call_nichede_function(
        "Filter_NDE", object_path, output_path, args=args, kwargs=kwargs, rscript=rscript
    )


def merge_objects(
    output_path: str | Path,
    *,
    args: list[Any] | None = None,
    kwargs: dict[str, Any] | None = None,
    rscript: str = "Rscript",
) -> Path | None:
    """Wrap the upstream ``MergeObjects`` helper."""
    return call_nichede_function(
        "MergeObjects", None, output_path, args=args, kwargs=kwargs, rscript=rscript
    )


def niche_lr_cell(
    object_path: str | Path,
    output_path: str | Path,
    *,
    args: list[Any] | None = None,
    kwargs: dict[str, Any] | None = None,
    rscript: str = "Rscript",
) -> Path | None:
    """Wrap the upstream ``niche_LR_cell`` helper."""
    return call_nichede_function(
        "niche_LR_cell", object_path, output_path, args=args, kwargs=kwargs, rscript=rscript
    )


def niche_lr_spot(
    object_path: str | Path,
    output_path: str | Path,
    *,
    args: list[Any] | None = None,
    kwargs: dict[str, Any] | None = None,
    rscript: str = "Rscript",
) -> Path | None:
    """Wrap the upstream ``niche_LR_spot`` helper."""
    return call_nichede_function(
        "niche_LR_spot", object_path, output_path, args=args, kwargs=kwargs, rscript=rscript
    )


def _generic_nichede_helper(
    function_name: str,
    output_path: str | Path,
    *,
    object_path: str | Path | None = None,
    args: list[Any] | None = None,
    kwargs: dict[str, Any] | None = None,
    rscript: str = "Rscript",
) -> Path | None:
    """Call a named upstream helper with consistent Python naming."""
    return call_nichede_function(
        function_name,
        object_path,
        output_path,
        args=args,
        kwargs=kwargs,
        rscript=rscript,
    )


def celltype_level(output_path: str | Path, **kwargs: Any) -> Path | None:
    """Wrap the upstream ``celltype_level`` helper."""
    return _generic_nichede_helper("celltype_level", output_path, **kwargs)


def celltype_level_fisher(output_path: str | Path, **kwargs: Any) -> Path | None:
    """Wrap the upstream ``celltype_level_fisher`` helper."""
    return _generic_nichede_helper("celltype_level_fisher", output_path, **kwargs)


def gene_level(output_path: str | Path, **kwargs: Any) -> Path | None:
    """Wrap the upstream ``gene_level`` helper."""
    return _generic_nichede_helper("gene_level", output_path, **kwargs)


def gene_level_fisher(output_path: str | Path, **kwargs: Any) -> Path | None:
    """Wrap the upstream ``gene_level_fisher`` helper."""
    return _generic_nichede_helper("gene_level_fisher", output_path, **kwargs)


def check_colloc(output_path: str | Path, **kwargs: Any) -> Path | None:
    """Wrap the upstream ``check_colloc`` helper."""
    return _generic_nichede_helper("check_colloc", output_path, **kwargs)


def contrast_post(output_path: str | Path, **kwargs: Any) -> Path | None:
    """Wrap the upstream ``contrast_post`` helper."""
    return _generic_nichede_helper("contrast_post", output_path, **kwargs)


def nb_lik(output_path: str | Path, **kwargs: Any) -> Path | None:
    """Wrap the upstream ``nb_lik`` helper."""
    return _generic_nichede_helper("nb_lik", output_path, **kwargs)


def t_to_p(output_path: str | Path, **kwargs: Any) -> Path | None:
    """Wrap the upstream ``T_to_p`` helper."""
    return _generic_nichede_helper("T_to_p", output_path, **kwargs)


def ultosymmetric(output_path: str | Path, **kwargs: Any) -> Path | None:
    """Wrap the upstream ``ultosymmetric`` helper."""
    return _generic_nichede_helper("ultosymmetric", output_path, **kwargs)


def call_nichede_function(
    function_name: str,
    object_path: str | Path | None = None,
    output_path: str | Path | None = None,
    *,
    args: list[Any] | None = None,
    kwargs: dict[str, Any] | None = None,
    rscript: str = "Rscript",
) -> Path | None:
    """Call an arbitrary exported NicheDE R function.

    This is an escape hatch for NicheDE helper functions not exposed as a
    first-class Python wrapper. When ``object_path`` is supplied, it is loaded as
    the first positional argument. If ``output_path`` is supplied, the R return
    value is saved as CSV for data frames and RDS otherwise.

    Args:
        function_name: Name of the R function in the ``nicheDE`` namespace.
        object_path: Optional RDS object passed as the first argument.
        output_path: Optional output path.
        args: JSON-serializable positional arguments.
        kwargs: JSON-serializable keyword arguments.
        rscript: Rscript executable.

    Returns:
        Output path when provided, otherwise ``None``.
    """
    output = Path(output_path).expanduser().resolve() if output_path else None
    run_r_bridge(
        "call_function",
        {
            "function_name": function_name,
            "object_path": str(Path(object_path).expanduser().resolve())
            if object_path is not None
            else None,
            "output_path": str(output) if output is not None else None,
            "args": args or [],
            "kwargs": kwargs or {},
        },
        rscript=rscript,
    )
    return output


def run_pipeline(
    adata: ad.AnnData,
    *,
    labels_key: str,
    spatial_key: str = "spatial",
    sigma: list[float] | None = None,
    layer: str | None = None,
    output_dir: str | Path = "outputs/nichede",
    num_cores: int = 1,
    large_scale: bool = False,
    batch_size: int = 1000,
    cutoff: float = 0.05,
    rscript: str = "Rscript",
    niche_de_kwargs: dict[str, Any] | None = None,
) -> NicheDEResult:
    """Run the standard NicheDE workflow from AnnData.

    Args:
        adata: Input AnnData object.
        labels_key: Column in ``adata.obs`` containing cell-type labels.
        spatial_key: Key in ``adata.obsm`` containing spatial coordinates.
        sigma: Kernel bandwidths. Defaults to ``[20, 200, 500]``.
        layer: Optional layer containing counts.
        output_dir: Output directory for inputs and results.
        num_cores: Number of cores for ``niche_DE``.
        large_scale: Use large-scale effective niche calculation.
        batch_size: Batch size for large-scale effective niche calculation.
        cutoff: Cutoff for large-scale effective niche calculation.
        rscript: Rscript executable.
        niche_de_kwargs: Additional keyword arguments for ``run_niche_de``.

    Returns:
        Paths produced by the workflow.
    """
    config = NicheDEConfig(
        labels_key=labels_key,
        spatial_key=spatial_key,
        sigma=sigma or [20, 200, 500],
        layer=layer,
        output_dir=output_dir,
        num_cores=num_cores,
        rscript=rscript,
        large_scale=large_scale,
        batch_size=batch_size,
        cutoff=cutoff,
    )

    input_dir = write_nichede_inputs(
        adata,
        config.resolved_output_dir() / "inputs",
        labels_key=config.labels_key,
        spatial_key=config.spatial_key,
        layer=config.layer,
    )
    library_matrix_path = create_library_matrix(input_dir, rscript=config.rscript)
    object_path = create_niche_de_object(
        input_dir,
        config.sigma,
        library_matrix_path=library_matrix_path,
        rscript=config.rscript,
    )
    effective_path = calculate_effective_niche(
        object_path,
        large_scale=config.large_scale,
        batch_size=config.batch_size,
        cutoff=config.cutoff,
        rscript=config.rscript,
    )
    results_path = run_niche_de(
        effective_path,
        num_cores=config.num_cores,
        rscript=config.rscript,
        **(niche_de_kwargs or {}),
    )
    return NicheDEResult(
        input_dir=input_dir,
        library_matrix_path=library_matrix_path,
        object_path=object_path,
        effective_niche_object_path=effective_path,
        niche_de_object_path=results_path,
    )
