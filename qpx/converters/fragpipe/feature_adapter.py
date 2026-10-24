"""FragPipe Feature adapter -- combined_ion.tsv to feature.parquet.

Loads a FragPipe ``combined_ion.tsv`` into DuckDB,
transforms feature rows into ``FeatureSchema``, and writes through ``FeatureWriter``.
"""

from __future__ import annotations

import logging
from typing import Optional

import pandas as pd

from qpx.converters.base import BaseConverter, resolve_columns
from qpx.converters.channel_labels import experiment_runs_from_sdrf
from qpx.converters.fragpipe.constants import is_decoy_accession, to_modifications, to_proforma
from qpx.converters.fragpipe.pg_adapter import (
    load_experiment_annotation,
    normalize_experiment_mapping,
    normalize_run_name,
)
from qpx.converters.mappings import get_field_mappings
from qpx.converters.utils import parse_uniprot_id, safe_float
from qpx.core.sql import escape_path, sql_build
from qpx.writers.feature import FeatureWriter

logger = logging.getLogger(__name__)

_PSM_LOOKUP_SQL = """
            SELECT
                CAST("Spectrum" AS VARCHAR) AS spectrum,
                CAST("Peptide" AS VARCHAR) AS sequence,
                COALESCE(CAST("Assigned Modifications" AS VARCHAR), '') AS assigned_mods,
                CAST("Charge" AS INTEGER) AS charge,
                $obs_expr AS obs_mz,
                $calc_expr AS calc_mz,
                $pp_expr AS pp_prob,
                CAST("Protein" AS VARCHAR) AS protein,
                $im_expr AS ion_mobility,
                $mc_expr AS missed_cleavages
            FROM _fp_psm_lookup"""

# Derive field map from central YAML mappings
_FEATURE_MAP = get_field_mappings("fragpipe", "feature")
_FEATURE_IDENTITY_COMPOSITE = (
    "quantification_unit_id",
    "peptidoform",
    "charge",
    "compensation_voltage",
)


def _extract_pg_proteins(protein_str: str) -> list[dict]:
    """Extract pg_protein structs from a FragPipe Protein field.

    Handles formats like ``sp|P12345|PROT_HUMAN`` or plain accession ``P12345``.
    Returns list of {accession, start, end, pre, post}; FragPipe gives no
    positions, so those are null.
    """
    if not protein_str:
        return []
    result = []
    for part in protein_str.split(","):
        part = part.strip()
        if not part:
            continue
        acc = parse_uniprot_id(part)[0]
        if acc:
            result.append(
                {
                    "accession": acc,
                    "start": None,
                    "end": None,
                    "pre": None,
                    "post": None,
                }
            )
    return result


class FragPipeFeatureAdapter(BaseConverter):
    """Convert FragPipe ``combined_ion.tsv`` to ``feature.parquet``.

    Usage::

        with FragPipeFeatureAdapter() as adapter:
            adapter.convert(
                feature_path="combined_ion.tsv",
                output_path="feature.parquet",
            )
    """

    def convert(
        self,
        feature_path: str,
        output_path: str,
        sdrf_path: Optional[str] = None,
        psm_path: Optional[str] = None,
        experiment_to_runs: dict[str, list[str]] | None = None,
        experiment_annotation_path: Optional[str] = None,
        chunksize: int = 500_000,
        creator: str = "fragpipe",
    ) -> None:
        """Run the combined_ion.tsv -> feature.parquet conversion.

        Args:
            feature_path: Path to FragPipe ``combined_ion.tsv``.
            output_path: Destination Parquet path.
            sdrf_path: Optional SDRF fallback for mapping experiments to raw runs.
            psm_path: Optional path to FragPipe ``psm.tsv`` for PSM-level lookups
                (mass error, PEP, scan, is_decoy).
            experiment_to_runs: Optional explicit mapping from each experiment
                intensity prefix to its member raw runs.
            experiment_annotation_path: Optional FragPipe
                ``experiment_annotation.tsv`` providing that mapping.
            chunksize: Rows per batch.
            creator: Creator tag in Parquet metadata.
        """
        # Step 1: Load feature file into DuckDB
        self._load_feature_file(feature_path)

        # Step 2: Resolve column mappings against actual input columns
        actual_cols = {
            c[0]
            for c in self._conn.execute(
                "SELECT column_name FROM information_schema.columns WHERE table_name='fragpipe_features'"
            ).fetchall()
        }
        self._resolved = resolve_columns(_FEATURE_MAP, actual_cols)

        # Step 3: Detect format and experiment columns
        self._require_ion_format()
        experiments = self._detect_experiment_columns()
        self.logger.info("Detected combined_ion format, experiments: %s", experiments)

        # Step 4: Build PSM lookup if psm.tsv provided
        psm_lookup = self._build_psm_lookup(psm_path) if psm_path else {}
        experiment_psm_sources = self._resolve_experiment_psm_sources(
            experiments,
            experiment_to_runs,
            experiment_annotation_path,
            sdrf_path,
            warn_on_ambiguous=bool(psm_lookup),
        )

        # Step 5: Stream and transform
        with FeatureWriter(
            output_path,
            creator=creator,
            compression=self._compression,
            identity_composite=_FEATURE_IDENTITY_COMPOSITE,
        ) as writer:
            for batch in self._query_batched("SELECT * FROM fragpipe_features", chunksize):
                df = batch.to_pandas()
                records = self._transform_batch(df, experiments, psm_lookup, experiment_psm_sources)
                if records:
                    self._track_scores(records)
                    writer.write_batch(records)

        self.logger.info(f"FragPipe feature conversion complete -> {output_path}")

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _load_feature_file(self, path: str) -> None:
        """Load combined_ion.tsv into DuckDB."""
        self._conn.execute(
            sql_build(
                """CREATE TABLE fragpipe_features AS
            SELECT * FROM read_csv_auto('$path',
                delim='\t', header=true, auto_detect=true)""",
                path=escape_path(path),
            )
        )
        count = self._conn.execute("SELECT COUNT(*) FROM fragpipe_features").fetchall()[0][0]
        self.logger.info(f"Loaded {count:,} FragPipe feature rows")

    def _require_ion_format(self) -> None:
        """Reject peptide-level input whose quantities cannot be split by charge."""
        cols = self._conn.execute(
            "SELECT column_name FROM information_schema.columns WHERE table_name='fragpipe_features'"
        ).fetchall()
        col_names = {c[0] for c in cols}
        if "Charge" in col_names:
            return
        if "Charges" in col_names:
            raise ValueError(
                "FragPipe combined_peptide.tsv is peptide-level: its intensity cannot be assigned to individual "
                "charges without duplicating quantity. Use combined_ion.tsv for Feature conversion."
            )
        raise ValueError("FragPipe Feature input must be combined_ion.tsv with a Charge column")

    def _detect_experiment_columns(self) -> list[str]:
        """Detect per-experiment intensity columns.

        FragPipe feature files have columns like ``<experiment> Intensity``.
        ``<experiment> MaxLFQ Intensity`` is a derived value for the same
        experiment, not a separate experiment.
        """
        cols = self._conn.execute(
            "SELECT column_name FROM information_schema.columns WHERE table_name='fragpipe_features'"
        ).fetchall()
        col_names = [c[0] for c in cols]

        experiments = set()
        for col in col_names:
            if col.endswith(" Intensity") and not col.endswith(" MaxLFQ Intensity"):
                experiments.add(col[: -len(" Intensity")])

        return sorted(experiments)

    def _resolve_experiment_psm_sources(
        self,
        experiments: list[str],
        experiment_to_runs: dict[str, list[str]] | None,
        experiment_annotation_path: str | None,
        sdrf_path: str | None,
        *,
        warn_on_ambiguous: bool,
    ) -> dict[str, str]:
        """Resolve experiments that can be enriched from exactly one raw run."""
        if experiment_to_runs is not None and experiment_annotation_path is not None:
            raise ValueError("Provide either experiment_to_runs or experiment_annotation_path, not both")
        if experiment_annotation_path is not None:
            experiment_runs = load_experiment_annotation(experiment_annotation_path)
        else:
            experiment_runs = normalize_experiment_mapping(experiment_to_runs or {})
        if not experiment_runs:
            experiment_runs = experiment_runs_from_sdrf(sdrf_path) or {}

        sources: dict[str, str] = {}
        ambiguous: list[str] = []
        for experiment in experiments:
            runs = experiment_runs.get(experiment)
            if runs is None:
                sources[experiment] = normalize_run_name(experiment)
            elif len(runs) == 1:
                sources[experiment] = normalize_run_name(runs[0])
            else:
                ambiguous.append(experiment)
        if warn_on_ambiguous and ambiguous:
            self.logger.warning(
                "Skipping PSM enrichment for multi-run FragPipe experiment(s): %s",
                ", ".join(ambiguous),
            )
        return sources

    # ------------------------------------------------------------------
    # PSM lookup
    # ------------------------------------------------------------------

    def _build_psm_lookup(self, psm_path: str) -> dict[tuple, dict]:
        """Build a lookup for feature-level PSM enrichment.

        Loads FragPipe ``psm.tsv`` into a temporary DuckDB table and extracts the
        first-occurrence PSM providing mass error, PEP, scan, and decoy flag.

        Each PSM is registered under its exact
        ``(source_file, peptidoform, charge)`` key. The experiment-to-source
        mapping is resolved separately from FragPipe metadata; no run-agnostic
        precursor fallback is used because it can attach one run's PSM to every
        experiment in a combined results table.
        """
        lookup: dict[tuple, dict] = {}

        self._conn.execute(
            sql_build(
                """CREATE TEMPORARY TABLE _fp_psm_lookup AS
            SELECT * FROM read_csv_auto('$path',
                delim='\t', header=true, auto_detect=true)""",
                path=escape_path(psm_path),
            )
        )
        count = self._conn.execute("SELECT COUNT(*) FROM _fp_psm_lookup").fetchall()[0][0]
        self.logger.info(f"PSM lookup: loaded {count:,} PSM rows from {psm_path}")

        # Detect columns
        actual_cols = {
            c[0]
            for c in self._conn.execute(
                "SELECT column_name FROM information_schema.columns WHERE table_name='_fp_psm_lookup'"
            ).fetchall()
        }

        sql = self._psm_lookup_query(actual_cols)
        rows = self._conn.execute(sql).fetchall()

        for values in rows:
            self._add_psm_lookup_row(lookup, values)

        # Clean up temporary table
        self._conn.execute("DROP TABLE IF EXISTS _fp_psm_lookup")
        self.logger.info(f"PSM lookup: {len(lookup):,} unique keys built")

        return lookup

    # ------------------------------------------------------------------
    # Transform
    # ------------------------------------------------------------------

    def _transform_batch(
        self,
        df: pd.DataFrame,
        experiments: list[str],
        psm_lookup: dict[tuple, dict] | None = None,
        experiment_psm_sources: dict[str, str] | None = None,
    ) -> list[dict]:
        # Pre-extract column arrays for faster per-row access than to_dict("records")
        col_arrays = {col: df[col].values for col in df.columns}
        rows = ({col: vals[i] for col, vals in col_arrays.items()} for i in range(len(df)))
        return self._transform_rows_with_skip(
            rows,
            lambda row: self._transform_row(row, experiments, psm_lookup, experiment_psm_sources),
            "FragPipe feature",
        )

    def _transform_row(
        self,
        row,
        experiments: list[str],
        psm_lookup: dict[tuple, dict] | None = None,
        experiment_psm_sources: dict[str, str] | None = None,
    ) -> list[dict]:
        """Transform a single row into one or more feature records.

        Expands into one record per experiment with non-zero intensity.
        """
        records: list[dict] = []
        common_fields, mz, protein_raw = self._feature_identity_fields(row)

        # Expand per experiment
        for experiment in experiments:
            int_col = f"{experiment} Intensity"
            intensity_val = safe_float(row.get(int_col)) or 0.0
            if intensity_val <= 0:
                continue

            intensities = [{"label": "LFQ", "intensity": float(intensity_val)}]

            # PSM lookup: enrich only through the exact raw source resolved for
            # this experiment. Multi-run experiments intentionally have no
            # source here because one combined Feature cannot represent one
            # arbitrarily selected fraction's PSM.
            psm_info = {}
            psm_source = (experiment_psm_sources or {}).get(experiment)
            if psm_lookup and psm_source:
                psm_info = psm_lookup.get((psm_source, common_fields["peptidoform"], str(common_fields["charge"])), {})

            psm_fields = self._feature_psm_fields(psm_info, mz, protein_raw)

            rec = {
                **common_fields,
                **psm_fields,
                "additional_scores": None,
                "predicted_rt": None,
                "run_file_name": experiment,
                "quantification_unit_id": experiment,
                "cv_params": None,
                "rt": None,
                "intensities": intensities,
                "additional_intensities": None,
                "pg_global_qvalue": None,
                "ion_mobility_start": None,
                "ion_mobility_stop": None,
                "gg_accessions": None,
                "id_run_file_name": experiment,
                "rt_start": None,
                "rt_stop": None,
            }
            records.append(rec)

        return records

    @staticmethod
    def _psm_lookup_query(actual_cols: set[str]) -> str:
        """Build the typed PSM enrichment query."""
        # Observed m/z: prefer calibrated
        obs_col = "Calibrated Observed M/Z" if "Calibrated Observed M/Z" in actual_cols else "Observed M/Z"
        calc_col = "Calculated M/Z"
        has_obs = obs_col in actual_cols
        has_calc = calc_col in actual_cols

        obs_expr = f'TRY_CAST("{obs_col}" AS DOUBLE)' if has_obs else "NULL"
        calc_expr = f'TRY_CAST("{calc_col}" AS DOUBLE)' if has_calc else "NULL"

        # PeptideProphet Probability → PEP = 1 - prob
        has_pp = "PeptideProphet Probability" in actual_cols
        pp_expr = 'TRY_CAST("PeptideProphet Probability" AS DOUBLE)' if has_pp else "NULL"

        # Ion mobility
        has_im = "Ion Mobility" in actual_cols
        im_expr = 'TRY_CAST("Ion Mobility" AS DOUBLE)' if has_im else "NULL"

        # Missed cleavages
        has_mc = "Number of Missed Cleavages" in actual_cols
        mc_expr = 'TRY_CAST("Number of Missed Cleavages" AS INTEGER)' if has_mc else "NULL"

        sql = sql_build(
            _PSM_LOOKUP_SQL,
            obs_expr=obs_expr,
            calc_expr=calc_expr,
            pp_expr=pp_expr,
            im_expr=im_expr,
            mc_expr=mc_expr,
        )
        return sql

    @staticmethod
    def _psm_source_scan(spectrum) -> tuple[str, int]:
        """Run name and scan number from a FragPipe ``Spectrum`` value; the scan is 0 when it does not parse."""
        # Parse source file and scan from Spectrum column
        tokens = str(spectrum).split(".") if spectrum else []
        source_file = normalize_run_name(tokens[0]) if tokens else ""
        scan_number = 0
        if len(tokens) >= 2:
            try:
                scan_number = int(tokens[1].lstrip("0") or "0")
            except ValueError:
                pass
        return source_file, scan_number

    def _add_psm_lookup_row(self, lookup: dict, values: tuple) -> None:
        """Register the first PSM for an exact source and precursor key."""
        (
            spectrum,
            sequence,
            assigned_mods,
            charge,
            obs_mz,
            calc_mz,
            pp_prob,
            protein,
            ion_mobility,
            missed_cleavages_val,
        ) = values
        source_file, scan_number = self._psm_source_scan(spectrum)

        sequence = str(sequence) if sequence else ""
        assigned_mods = str(assigned_mods) if assigned_mods else ""
        peptidoform = to_proforma(assigned_mods, sequence)
        charge_str = str(charge) if charge is not None else "0"

        key = (source_file, peptidoform, charge_str)
        if key not in lookup:
            # PEP = 1 - PeptideProphet Probability
            pep = None
            if pp_prob is not None:
                pep = 1.0 - pp_prob if pp_prob <= 1.0 else pp_prob

            is_decoy = is_decoy_accession(protein)

            info = {
                "observed_mz": float(obs_mz) if obs_mz is not None else None,
                "calculated_mz": (float(calc_mz) if calc_mz is not None else None),
                "pep": pep,
                "is_decoy": is_decoy,
                "scan": [scan_number] if scan_number > 0 else [],
                "ion_mobility": (float(ion_mobility) if ion_mobility is not None else None),
                "missed_cleavages": (int(missed_cleavages_val) if missed_cleavages_val is not None else None),
            }
            lookup[key] = info

    def _feature_identity_fields(self, row: dict) -> tuple:
        """Parse the feature peptide and protein identity fields."""
        r = self._resolved  # shorthand for resolved column mappings

        # Sequence
        sequence = str(row.get(r.get("sequence", "Peptide Sequence"), ""))

        # Peptidoform -- build ProForma from sequence + Assigned Modifications
        mods_col = r.get("modifications", "Assigned Modifications")
        assigned_mods_raw = row.get(mods_col, "")
        assigned_mods_str = str(assigned_mods_raw) if pd.notna(assigned_mods_raw) and assigned_mods_raw else ""
        peptidoform = to_proforma(assigned_mods_str, sequence)

        # Protein mapping (schema: list<pg_protein>; FragPipe does not provide start/end)
        protein_raw = str(row.get(r.get("pg_accessions", "Protein"), ""))
        pg_accessions = _extract_pg_proteins(protein_raw) or None
        anchor_protein = pg_accessions[0]["accession"] if pg_accessions else ""

        # Gene
        gene_raw = row.get(r.get("gg_names", "Gene"), "")
        gg_names = [str(gene_raw)] if pd.notna(gene_raw) and gene_raw else None

        # Modifications (reuse assigned_mods_str already extracted for peptidoform)
        modifications = None
        if assigned_mods_str:
            _, modifications = to_modifications(assigned_mods_str, sequence)

        # M/Z (from feature file — used as fallback)
        mz = safe_float(row.get(r.get("observed_mz", "M/Z")))
        compensation_voltage = safe_float(row.get(r.get("compensation_voltage", "Compensation Voltage"))) or 0.0

        charge = int(row.get(r.get("charge", "Charge"), 0))
        return (
            {
                "sequence": sequence,
                "peptidoform": peptidoform,
                "modifications": modifications,
                "charge": charge,
                "compensation_voltage": compensation_voltage,
                "pg_accessions": pg_accessions,
                "anchor_protein": anchor_protein,
                "unique": len(pg_accessions) <= 1 if pg_accessions else True,
                "gg_names": gg_names,
            },
            mz,
            protein_raw,
        )

    @staticmethod
    def _feature_psm_fields(psm_info: dict, mz: float | None, protein_raw: str) -> dict:
        """Resolve spectrum fields from exact PSM enrichment or feature values."""
        _calc = psm_info.get("calculated_mz")
        _obs = psm_info.get("observed_mz")
        mass_error_ppm = 1e6 * (_obs - _calc) / _calc if _calc and _obs else None

        # Is decoy: prefer PSM lookup; fall back to the decoy prefix on the RAW
        # protein field (before parse_uniprot_id strips e.g. the "rev_" prefix).
        _is_decoy_fallback = is_decoy_accession(protein_raw)
        return {
            "posterior_error_probability": psm_info.get("pep"),
            "is_decoy": psm_info.get("is_decoy", _is_decoy_fallback),
            "calculated_mz": _calc or mz,
            "observed_mz": _obs or mz,
            "mass_error_ppm": mass_error_ppm,
            "scan": psm_info.get("scan", []),
            "ion_mobility": psm_info.get("ion_mobility"),
            "missed_cleavages": psm_info.get("missed_cleavages"),
        }
