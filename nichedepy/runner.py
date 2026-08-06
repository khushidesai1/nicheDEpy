"""Subprocess runner for the NicheDE R bridge."""

from __future__ import annotations

import json
import subprocess
import tempfile
from importlib.resources import files
from pathlib import Path
from typing import Any


def r_bridge_path() -> Path:
    """Return the bundled R bridge path."""
    return Path(files("nichedepy").joinpath("r/nichede_bridge.R"))


def run_r_bridge(
    command: str,
    params: dict[str, Any],
    *,
    rscript: str = "Rscript",
    check: bool = True,
) -> subprocess.CompletedProcess[str]:
    """Run a NicheDE bridge command through ``Rscript``.

    Args:
        command: Bridge command name.
        params: JSON-serializable command parameters.
        rscript: Rscript executable.
        check: Whether to raise an error when R exits with a non-zero status.

    Returns:
        The completed subprocess result.
    """
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as handle:
        json.dump(params, handle)
        params_path = Path(handle.name)

    try:
        cmd = [rscript, str(r_bridge_path()), command, str(params_path)]
        result = subprocess.run(
            cmd,
            check=False,
            text=True,
            capture_output=True,
        )
    finally:
        params_path.unlink(missing_ok=True)

    if check and result.returncode != 0:
        message = (
            "NicheDE R bridge failed.\n"
            f"Command: {command}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )
        raise RuntimeError(message)
    return result

