"""Transform CLI dataset implementation."""

from pathlib import Path
from typing import Optional

import click


def _copy_dataset_files(dataset: Path, out_dir: Path) -> None:
    """Copy the complete dataset, including directory-backed views."""
    import shutil

    if out_dir.exists() and any(out_dir.iterdir()):
        raise click.ClickException("--output-folder must be empty")
    shutil.copytree(dataset, out_dir, dirs_exist_ok=True)


def _gene_map_dataset_prefix(dataset: Path, command: str = "gene-map", prefix: Optional[str] = None) -> str:
    """Reject unsupported quantification layouts before copying or annotation.

    ``prefix`` skips discovery when the caller already knows it (a converter that
    just wrote the files), so a second dataset in the same folder cannot make the
    prefix ambiguous.
    """
    from qpx.transforms.utils import discover_qpx_file_prefix

    for view in ("pg", "feature"):
        if any((dataset / view).rglob("*.parquet")):
            raise click.ClickException(f"Partitioned {view} data is not supported by {command}; use flat Parquet views")
    if prefix:
        return prefix
    try:
        return discover_qpx_file_prefix(dataset)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc


def _stage_dataset_integrity(qpx_dataset, prefix: str, staging: Path, out_dir: Path) -> tuple[Path, Path] | None:
    """Update existing integrity records for the staged quantification files."""
    from qpx.dataset import Dataset
    from qpx.transforms.utils import count_staged_structures
    from qpx.writers import DatasetWriter

    if qpx_dataset.dataset_meta is None:
        return None
    metadata = qpx_dataset.dataset_meta.to_df()
    if metadata.empty:
        return None
    record = metadata.iloc[0].to_dict()
    fields = ("file_checksums", "file_row_counts", "file_sizes_bytes")
    if not any(isinstance(record.get(field), dict) for field in fields):
        return None

    with Dataset(str(staging), file_prefix=prefix) as staged_dataset:
        integrity = staged_dataset.compute_integrity()
    for field in fields:
        if isinstance(record.get(field), dict):
            record[field].update(integrity[field])
    record["total_structures"] = count_staged_structures(out_dir, staging, integrity["file_checksums"])
    record["packaged_at"] = integrity["packaged_at"]
    name = f"{prefix}.dataset.parquet"
    with DatasetWriter(staging / name) as writer:
        writer.write_batch([record])
    return staging / name, out_dir / name


def _refresh_dataset_mudata(out_dir: Path, prefix: str) -> None:
    """Rebuild the dataset's MuData view, so a stale h5mu never outlives the parquet."""
    from qpx.mudata import write_dataset_mudata

    if not (out_dir / f"{prefix}.h5mu").is_file():
        return

    refreshed = write_dataset_mudata(out_dir, prefix)
    if refreshed is None:
        click.echo("  h5mu: could not be rebuilt (see log); parquet views are annotated")
    else:
        click.echo(f"  {refreshed.name}: MuData view rebuilt")
