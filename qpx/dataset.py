"""Dataset class — the central entry point for opening QPX datasets."""

from __future__ import annotations

import logging
from functools import cached_property
from pathlib import Path

import duckdb

from qpx._version import __version__
from qpx.core.data import (
    PG,
    PSM,
    BaseStructure,
    DatasetMeta,
    Feature,
    MzSpectra,
    Ontology,
    PepMap,
    Proteoform,
    Provenance,
    Run,
    Sample,
    Site,
)
from qpx.core.dataset_queries import _DatasetQueryMixin
from qpx.core.dataset_validation import _DatasetValidationMixin
from qpx.core.engine import DuckDBEngine
from qpx.core.sql import sql_build, validate_table
from qpx.version import (
    QPX_SPEC_VERSION,
    QpxVersionError,
    check_feature_file_compatible,
    check_feature_schema_compatible,
    check_pg_columns_compatible,
    check_pg_file_compatible,
)

_log = logging.getLogger(__name__)

# Per-file layout guards run on every local file of a structure before it is registered.
_FILE_LAYOUT_GUARDS = {"pg": check_pg_file_compatible, "feature": check_feature_file_compatible}


class Dataset(_DatasetValidationMixin, _DatasetQueryMixin):
    """
    A QPX dataset — a directory of related Parquet/H5AD files.

    Usage:
        ds = qpx.Dataset("PXD014414/")
        ds.feature.filter("charge > 2").to_df()
        ds.pg.join(ds.run).to_df()

    Partial access:
        ds = qpx.Dataset("path/", structures=["feature", "sample", "run"])
    """

    # Data structure registry: name → (Class, file_suffix)
    _STRUCTURE_REGISTRY = {
        "psm": (PSM, ".psm.parquet"),
        "feature": (Feature, ".feature.parquet"),
        "pg": (PG, ".pg.parquet"),
        "mz": (MzSpectra, ".mz.parquet"),
        "sample": (Sample, ".sample.parquet"),
        "run": (Run, ".run.parquet"),
        "dataset": (DatasetMeta, ".dataset.parquet"),
        "ontology": (Ontology, ".ontology.parquet"),
        "provenance": (Provenance, ".provenance.parquet"),
        "pepmap": (PepMap, ".pepmap.parquet"),
        "site": (Site, ".site.parquet"),
        "proteoform": (Proteoform, ".proteoform.parquet"),
    }

    @classmethod
    def structure_names(cls) -> list[str]:
        """Every structure name QPX knows about, in registry order.

        The single source of truth for anything that needs to enumerate
        structures — CLI choices in particular. Hand-maintained copies drifted:
        ``pepmap`` was registered here and exposed on the API while both CLIs
        rejected it (bigbio/qpx#289).
        """
        return list(cls._STRUCTURE_REGISTRY)

    def __init__(
        self,
        path: str | Path,
        structures: list[str] | None = None,
        duckdb_memory: str | None = None,
        duckdb_threads: int | None = None,
        s3_config: dict | None = None,
        file_prefix: str | None = None,
    ):
        self._is_s3 = isinstance(path, str) and path.startswith("s3://")
        self.path = path if self._is_s3 else Path(path)
        self._file_prefix = file_prefix
        self._engine = DuckDBEngine(
            memory_limit=duckdb_memory,
            threads=duckdb_threads,
            s3_config=s3_config if self._is_s3 else None,
        )

        self._structures: dict[str, BaseStructure] = {}
        self._discover_and_register(structures)

    def _discover_and_register(self, requested: list[str] | None):
        """Scan directory, find QPX files, register as DuckDB tables.

        Checks for single Parquet files first, then falls back to
        Hive-partitioned directories (e.g., feature/ with run_file_name= subdirs).
        For S3 paths, attempts to register each structure via S3 glob.
        """
        if self._is_s3:
            self._discover_s3(requested)
        else:
            self._discover_local(requested)

    def _discover_s3(self, requested: list[str] | None):
        """Register structures from S3 path."""
        for name, (cls, suffix) in self._STRUCTURE_REGISTRY.items():
            if requested and name not in requested:
                continue
            file_name = f"{self._file_prefix}{suffix}" if self._file_prefix else f"*{suffix}"
            s3_glob = f"{self.path.rstrip('/')}/{file_name}"
            try:
                # PyArrow cannot open an S3 glob to read the footer, so register
                # first and run the PG guard against the DuckDB view's columns.
                self._engine.register_s3_parquet(name, s3_glob)
                if name == "pg":
                    cols = [row[0] for row in self._engine.execute('DESCRIBE "pg"').fetchall()]
                    check_pg_columns_compatible(cols, source=s3_glob)
                elif name == "feature":
                    schema = self._engine.connection.sql('SELECT * FROM "feature" LIMIT 0').to_arrow_table().schema
                    check_feature_schema_compatible(schema, s3_glob)
                self._structures[name] = cls(
                    engine=self._engine,
                    table_name=name,
                    file_path=f"{self.path}/{name}",
                    # The glob is what was actually registered; a consumer
                    # re-registering from the display path would read nothing.
                    file_paths=[s3_glob],
                )
            except QpxVersionError as exc:
                # One incompatible/old structure file must not abort the whole
                # dataset — skip it with a warning so the other structures stay
                # usable. (A direct read_pg()/PG.from_file() still raises.)
                self._engine.execute(
                    sql_build(
                        "DROP VIEW IF EXISTS $view",
                        view=validate_table(name),
                    )
                )
                _log.warning("Skipping incompatible structure '%s': %s", name, exc)
            except (FileNotFoundError, duckdb.IOException):
                pass  # Structure not present in S3
            except (OSError, ValueError, duckdb.Error) as exc:
                _log.warning("Failed to register '%s' from S3: %s", name, exc)

    def _discover_local(self, requested: list[str] | None):
        """Register structures from local filesystem.

        Each structure is registered independently: a bad/old/garbage file for
        one structure is skipped with a warning so the rest of the dataset stays
        usable, instead of aborting the whole Dataset construction.
        """
        for name, (cls, suffix) in self._STRUCTURE_REGISTRY.items():
            if requested and name not in requested:
                continue
            try:
                self._register_local_structure(name, cls, suffix)
            except QpxVersionError as exc:
                # One incompatible/old structure file must not abort the whole
                # dataset — skip it with a warning so the other structures stay
                # usable. (A direct read_pg()/PG.from_file() still raises.)
                _log.warning("Skipping incompatible structure '%s': %s", name, exc)
            except (OSError, ValueError, duckdb.Error) as exc:
                _log.warning("Failed to register structure '%s': %s", name, exc)

    def _register_local_structure(self, name: str, cls, suffix: str):
        """Register a single local structure (single file or Hive-partitioned dir)."""
        # Check for single file first
        pattern = f"{self._file_prefix}{suffix}" if self._file_prefix else f"*{suffix}"
        matches = sorted(self.path.glob(pattern))
        if matches:
            # Read ALL matching shards unioned, matching the S3 path (which globs
            # and reads every match). Previously only the first was taken, so a
            # sharded structure returned different rows locally vs over S3
            # (bigbio/qpx#252). The registered structure's file_path points at the
            # first shard for provenance; the view unions them all.
            guard = _FILE_LAYOUT_GUARDS.get(name)
            if guard is not None:
                for match in matches:
                    guard(match)
            self._engine.register_parquet_files(name, matches)
            self._structures[name] = cls(
                engine=self._engine,
                table_name=name,
                file_path=matches[0],
                file_paths=matches,
            )
        elif self._file_prefix is None:
            # Check for Hive-partitioned directory
            part_dir = self.path / name
            part_files = list(part_dir.glob("**/*.parquet")) if part_dir.is_dir() else []
            if part_files:
                guard = _FILE_LAYOUT_GUARDS.get(name)
                if guard is not None:
                    for part_file in sorted(part_files):
                        guard(part_file)
                self._engine.register_partitioned_parquet(name, part_dir)
                self._structures[name] = cls(
                    engine=self._engine,
                    table_name=name,
                    file_path=part_dir,
                )

    # --- Data structure accessors (lazy, return None if not present) ---
    @property
    def psm(self) -> PSM | None:
        return self._structures.get("psm")

    @property
    def feature(self) -> Feature | None:
        return self._structures.get("feature")

    @property
    def pg(self) -> PG | None:
        return self._structures.get("pg")

    @property
    def mz(self) -> MzSpectra | None:
        return self._structures.get("mz")

    @property
    def sample(self) -> Sample | None:
        return self._structures.get("sample")

    @property
    def run(self) -> Run | None:
        return self._structures.get("run")

    @property
    def dataset_meta(self) -> DatasetMeta | None:
        return self._structures.get("dataset")

    @property
    def ontology(self) -> Ontology | None:
        return self._structures.get("ontology")

    @property
    def provenance(self) -> Provenance | None:
        return self._structures.get("provenance")

    @property
    def pepmap(self) -> PepMap | None:
        return self._structures.get("pepmap")

    @property
    def site(self) -> Site | None:
        """The site view, or None when the dataset has no site.parquet."""
        return self._structures.get("site")

    @property
    def proteoform(self) -> Proteoform | None:
        """The proteoform view, or None when the dataset has no proteoform.parquet."""
        return self._structures.get("proteoform")

    # --- View accessors (cached) ---
    @cached_property
    def protein_view(self):
        from qpx.views.api import ProteinView

        return ProteinView(self)

    @cached_property
    def peptide_view(self):
        from qpx.views.api import PeptideView

        return PeptideView(self)

    @cached_property
    def identification_summary(self):
        from qpx.views.api import IdentificationSummaryView

        return IdentificationSummaryView(self)

    @cached_property
    def run_summary(self):
        from qpx.views.api import RunSummaryView

        return RunSummaryView(self)

    @cached_property
    def modification_view(self):
        from qpx.views.api import ModificationView

        return ModificationView(self)

    @cached_property
    def qc_view(self):
        from qpx.views.api import QualityControlView

        return QualityControlView(self)

    @cached_property
    def sample_summary(self):
        from qpx.views.api import SampleSummaryView

        return SampleSummaryView(self)

    @cached_property
    def ae_view(self):
        from qpx.views.api import AbsoluteExpressionView

        return AbsoluteExpressionView(self)

    @property
    def available_structures(self) -> list[str]:
        return list(self._structures.keys())

    @property
    def file_prefix(self) -> str | None:
        """Prefix of the QPX file names this dataset was opened with, if any."""
        return self._file_prefix

    # --- PRIDE enrichment ---
    def enrich_from_pride(self, project_accession: str | None = None) -> dict:
        """Fetch PRIDE metadata and update dataset.parquet.

        Retrieves project title, description, and PubMed ID from the
        PRIDE REST API and writes them into the dataset metadata.

        Parameters
        ----------
        project_accession : str or None
            PRIDE/ProteomeXchange accession. If *None*, reads from
            the existing dataset.parquet.

        Returns
        -------
        dict
            The fetched PRIDE metadata dict.

        Raises
        ------
        ValueError
            If no accession is provided and none is found in dataset.parquet,
            or if the accession is not found in PRIDE.
        """
        self._require_local("enrich_from_pride")

        # Resolve accession
        if project_accession is None:
            if self.dataset_meta is None:
                raise ValueError("No project_accession provided and no dataset.parquet found.")
            meta_df = self.dataset_meta.to_df()
            project_accession = meta_df["project_accession"].iloc[0]
            if not project_accession or project_accession == "unknown":
                raise ValueError("No valid project_accession in dataset.parquet. Provide one explicitly.")

        from qpx.core.pride import fetch_pride_metadata

        metadata = fetch_pride_metadata(project_accession)

        # Read existing dataset record, update fields, and rewrite
        if self.dataset_meta is not None:
            meta_df = self.dataset_meta.to_df()
            record = meta_df.iloc[0].to_dict()
        else:
            record = {
                "project_accession": project_accession,
                "creation_date": None,
            }

        record["project_title"] = metadata["project_title"]
        record["project_description"] = metadata["project_description"]
        record["pubmed_id"] = metadata["pubmed_id"]

        # Derive prefix from existing dataset file to overwrite in-place
        prefix = None
        if self.dataset_meta is not None:
            ds_suffix = self._STRUCTURE_REGISTRY["dataset"][1]
            existing_name = Path(self.dataset_meta._file_path).name
            if existing_name.endswith(ds_suffix):
                prefix = existing_name[: -len(ds_suffix)]

        self.save_structure([record], "dataset", prefix=prefix)
        self.refresh()

        return metadata

    # --- Write-back ---
    # Writer registry: name → (WriterClassName, file_suffix)
    _WRITER_REGISTRY = {
        "psm": ("PsmWriter", ".psm.parquet"),
        "feature": ("FeatureWriter", ".feature.parquet"),
        "pg": ("PgWriter", ".pg.parquet"),
        "mz": ("MzWriter", ".mz.parquet"),
        "sample": ("SampleWriter", ".sample.parquet"),
        "run": ("RunWriter", ".run.parquet"),
        "dataset": ("DatasetWriter", ".dataset.parquet"),
        "ontology": ("OntologyWriter", ".ontology.parquet"),
        "provenance": ("ProvenanceWriter", ".provenance.parquet"),
        "pepmap": ("PepMapWriter", ".pepmap.parquet"),
        "site": ("SiteWriter", ".site.parquet"),
        "proteoform": ("ProteoformWriter", ".proteoform.parquet"),
    }

    def save_structure(
        self,
        data,
        structure: str,
        prefix: str | None = None,
    ) -> Path:
        """Write validated Parquet back into the dataset directory.

        Only supported for local datasets. For S3-backed datasets, use the
        appropriate writer (e.g. ``FeatureWriter(path)``) with an explicit
        local output path.

        Parameters
        ----------
        data : list[dict] | pd.DataFrame | pa.Table
            Records to write. Schema-validated by the QPX writer.
        structure : str
            Structure name (e.g., "feature", "pg", "sample").
        prefix : str | None
            File name prefix. Defaults to the dataset directory name.

        Returns
        -------
        Path
            Path to the written Parquet file.
        """
        import pyarrow as pa

        self._require_local("save_structure")
        if structure not in self._WRITER_REGISTRY:
            raise ValueError(f"Unknown structure '{structure}'. Valid: {list(self._WRITER_REGISTRY.keys())}")

        writer_name, suffix = self._WRITER_REGISTRY[structure]

        import qpx.writers as writers_mod

        writer_cls = getattr(writers_mod, writer_name)

        path = Path(self.path)
        file_prefix = prefix or path.name
        output_path = path / f"{file_prefix}{suffix}"

        with writer_cls(output_path) as writer:
            if isinstance(data, list):
                writer.write_batch(data)
            elif isinstance(data, pa.Table):
                # Project onto the writer's schema (adds any absent columns such
                # as a derived id or optional cross-refs, then casts to normalize
                # nullability/types). write_table stamps the derived id.
                writer.write_table(writer.align_table_to_schema(data))
            else:
                writer.write_dataframe(data)

        return output_path

    # AnnData view type -> canonical file suffix (before .h5ad)
    _ANNDATA_VIEWS = {"ae", "de"}

    def save_anndata(
        self,
        adata,
        name: str | None = None,
        *,
        view: str | None = None,
    ) -> Path:
        """Write an AnnData object to the dataset directory as ``.h5ad``.

        Only supported for local datasets. For S3-backed datasets, write to
        a local path with ``adata.write_h5ad(path)``.

        Parameters
        ----------
        adata : anndata.AnnData
            The AnnData object to save.
        name : str, optional
            Base name for the file (without extension), e.g., "de_results".
            When omitted the dataset prefix is used.
        view : str, optional
            AnnData view type: ``"ae"`` (absolute expression) or ``"de"``
            (differential expression).  When provided the output file
            follows the QPX naming convention
            ``<prefix>.<view>.h5ad`` (e.g., ``PXD000000.ae.h5ad``).
            Ignored when *name* is explicitly given.

        Returns
        -------
        Path
            Path to the written .h5ad file.
        """
        self._require_local("save_anndata")
        path = Path(self.path)
        if name is not None:
            output_path = path / f"{name}.h5ad"
        elif view is not None:
            if view not in self._ANNDATA_VIEWS:
                raise ValueError(f"Unknown AnnData view '{view}'. Choose from: {sorted(self._ANNDATA_VIEWS)}")
            prefix = path.name
            output_path = path / f"{prefix}.{view}.h5ad"
        else:
            raise ValueError("Either 'name' or 'view' must be provided")
        adata.uns["qpx_version"] = QPX_SPEC_VERSION
        adata.uns["writer_version"] = __version__
        adata.write_h5ad(output_path)
        return output_path

    # --- External registration (for downstream tool integration) ---
    def register_external(self, name: str, file_path: str | Path) -> None:
        """Register an external Parquet file in the DuckDB engine."""
        self._engine.register_parquet(name, file_path)

    # --- Lifecycle ---
    def refresh(self) -> None:
        """Re-scan the dataset directory for new or updated files.

        Clears all cached view instances so next access gets fresh data.
        Call this after writing new structures or AnnData files.
        """
        for attr in (
            "protein_view",
            "peptide_view",
            "identification_summary",
            "run_summary",
            "modification_view",
            "qc_view",
            "sample_summary",
            "ae_view",
        ):
            self.__dict__.pop(attr, None)

        self._structures.clear()
        self._discover_and_register(None)

    def close(self):
        self._engine.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def __repr__(self):
        structs = ", ".join(self.available_structures)
        return f"Dataset('{self.path}', structures=[{structs}])"
