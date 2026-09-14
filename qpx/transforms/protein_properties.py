"""Fill protein properties from a FASTA when the producer did not record them.

Some producers never report a protein's sequence (DIA-NN), and some consensusXML
files carry protein hits without it (MSV000085836). Given the FASTA used for the
search, this module fills, for target rows only and only where the value is
null:

- ``pg.molecular_weight``   — average mass of the anchor protein, in kDa
- ``pg.sequence_coverage``  — percent of the anchor covered by the dataset's
  peptides mapped to it
- ``feature.pg_positions``  — every one-based occurrence of the peptide in each
  member of its protein group

A value the producer recorded is never overwritten. The FASTA is optional: a
protein absent from it (a DIA-NN internal decoy, a contaminant from another
database, an isoform not in this FASTA) keeps a null value, and the result reports
how many rows could not be matched so a wrong FASTA is visible rather than
silently partial.
"""

from __future__ import annotations

import gzip
import hashlib
import logging
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from qpx.core.protein_sequence import (
    average_molecular_weight_kda,
    peptide_occurrences,
    sequence_coverage_percent,
)

logger = logging.getLogger(__name__)

# Union of the decoy conventions used across qpx converters and quantms databases.
_DECOY_PREFIXES = ("DECOY", "REV_", "RANDOM_", "XXX_")


def _is_decoy_identifier(identifier: str) -> bool:
    """True for a FASTA entry whose identifier, or its accession field, is a decoy."""
    upper = identifier.upper()
    if upper.startswith(_DECOY_PREFIXES):
        return True
    parts = identifier.split("|")
    return len(parts) >= 2 and parts[1].upper().startswith(_DECOY_PREFIXES)


class FastaSequences:
    """Protein sequences from a FASTA, looked up by full identifier or accession.

    ``>sp|P12345|NAME ...`` is reachable as ``sp|P12345|NAME`` (OpenMS) and as
    ``P12345`` (DIA-NN). Decoy entries are skipped: a quantms target-decoy database
    carries ``DECOY_sp|P12345|NAME``, whose accession field would otherwise collide
    with the target and erase it as a conflict. A key that maps to two different
    sequences is dropped rather than guessed.
    """

    def __init__(self) -> None:
        self._sequences: dict[str, str] = {}
        self._conflicting: set[str] = set()
        self.entries = 0
        self.decoy_entries = 0

    @classmethod
    def from_path(cls, path: str | Path) -> FastaSequences:
        """Index a FASTA (plain or ``.gz``); multi-line sequences are joined."""
        index = cls()
        opener = gzip.open if str(path).endswith(".gz") else open
        identifier: str | None = None
        chunks: list[str] = []
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith(">"):
                    index._add(identifier, chunks)
                    header = line[1:].strip()
                    identifier = header.split(None, 1)[0] if header else ""
                    chunks = []
                elif identifier is not None:
                    chunks.append(line.strip())
        index._add(identifier, chunks)
        if index._conflicting:
            logger.warning(
                "%d FASTA identifier(s) map to different sequences and are left unresolved, e.g. %s",
                len(index._conflicting),
                sorted(index._conflicting)[:3],
            )
        return index

    def _add(self, identifier: str | None, chunks: list[str]) -> None:
        if not identifier:
            return
        self.entries += 1
        if _is_decoy_identifier(identifier):
            self.decoy_entries += 1
            return
        sequence = "".join(chunks).upper().rstrip("*")
        if not sequence:
            return
        keys = {identifier}
        parts = identifier.split("|")
        if len(parts) >= 2 and parts[1]:
            keys.add(parts[1])
        for key in keys:
            if key in self._conflicting:
                continue
            known = self._sequences.get(key)
            if known is None:
                self._sequences[key] = sequence
            elif known != sequence:
                del self._sequences[key]
                self._conflicting.add(key)

    def __len__(self) -> int:
        return len(self._sequences)

    def get(self, accession: str | None) -> str | None:
        """The sequence for ``accession``, or None when absent or ambiguous."""
        if not accession or accession in self._conflicting:
            return None
        sequence = self._sequences.get(accession)
        if sequence is None and "|" in accession:
            parts = accession.split("|")
            if len(parts) >= 2 and parts[1] not in self._conflicting:
                sequence = self._sequences.get(parts[1])
        return sequence


@dataclass
class ProteinPropertiesReport:  # pylint: disable=too-many-instance-attributes
    """What was filled, and what could not be matched to the FASTA (plain counters)."""

    fasta_entries: int = 0
    fasta_decoy_entries: int = 0
    pg_rows: int = 0
    pg_rows_eligible: int = 0
    pg_coverage_filled: int = 0
    pg_molecular_weight_filled: int = 0
    pg_anchors_not_in_fasta: int = 0
    feature_rows: int = 0
    feature_rows_eligible: int = 0
    feature_positions_filled: int = 0
    feature_rows_without_match: int = 0
    unmatched_examples: list[str] = field(default_factory=list)

    @property
    def pg_match_rate(self) -> float | None:
        """Share of eligible target pg rows whose anchor was found in the FASTA."""
        if not self.pg_rows_eligible:
            return None
        return 1.0 - self.pg_anchors_not_in_fasta / self.pg_rows_eligible


def _note_unmatched(report: ProteinPropertiesReport, accession: str) -> None:
    if len(report.unmatched_examples) < 5 and accession not in report.unmatched_examples:
        report.unmatched_examples.append(accession)


def _accession_lists(column: pa.ChunkedArray) -> list[tuple[str, ...] | None]:
    """Per-row accession tuples from a ``list<struct<accession,...>>`` column.

    Vectorised: flattening the struct field avoids materialising every struct
    as a Python dict, which dominates the cost on a 20M-row feature view.
    """
    import pyarrow.compute as pc

    out: list[tuple[str, ...] | None] = []
    for chunk in column.chunks:
        if len(chunk) == 0:
            continue
        offsets = chunk.offsets.to_numpy()
        accessions = pc.struct_field(chunk.flatten(), "accession").to_pylist()  # pylint: disable=no-member
        validity = chunk.is_valid().to_pylist()
        base = offsets[0]
        for row, valid in enumerate(validity):
            if not valid:
                out.append(None)
                continue
            lo, hi = offsets[row] - base, offsets[row + 1] - base
            out.append(tuple(a for a in accessions[lo:hi] if a))
    return out


def _peptide_protein_candidates(feature_path: Path) -> dict[str, set[str]]:
    """Map each target peptide sequence to the proteins it is attributed to.

    From the feature view's group memberships, plus the proteins named by
    positions a producer already recorded (a peptide shared across groups has a
    null group but recorded positions).
    """
    import duckdb

    con = duckdb.connect()
    path = str(feature_path).replace("'", "''")
    columns = {row[0] for row in con.execute(f"DESCRIBE SELECT * FROM read_parquet('{path}')").fetchall()}
    decoy = "NOT coalesce(is_decoy, false)" if "is_decoy" in columns else "true"
    candidates: dict[str, set[str]] = defaultdict(set)
    queries = [
        f"""SELECT sequence, list_distinct(list_transform(pg_accessions, x -> x.accession))
            FROM read_parquet('{path}') WHERE {decoy} AND sequence IS NOT NULL AND pg_accessions IS NOT NULL
            GROUP BY ALL"""
    ]
    if "pg_positions" in columns:
        queries.append(
            f"""SELECT sequence, list_distinct(list_transform(pg_positions, x -> x.protein_accession))
                FROM read_parquet('{path}') WHERE {decoy} AND sequence IS NOT NULL AND pg_positions IS NOT NULL
                GROUP BY ALL"""
        )
    for query in queries:
        for sequence, accessions in con.execute(query).fetchall():
            for accession in accessions or ():
                if accession:
                    candidates[sequence].add(accession)
    con.close()
    return candidates


def _protein_peptides(candidates: dict[str, set[str]]) -> dict[str, set[str]]:
    proteins: dict[str, set[str]] = defaultdict(set)
    for sequence, accessions in candidates.items():
        for accession in accessions:
            proteins[accession].add(sequence)
    return proteins


def _stamped_schema(schema: pa.Schema) -> tuple[pa.Schema, str]:
    """The source schema with its footer re-stamped (new uuid/date, same identity)."""
    from qpx.writers.base import _stamp_footer_metadata

    metadata = schema.metadata or {}
    compression = metadata.get(b"compression_format", b"zstd").decode() or "zstd"
    stamped = _stamp_footer_metadata(schema.empty_table(), compression).schema
    return stamped, compression


def _rewrite_view(source: Path, destination: Path, fill_batch) -> None:
    """Rewrite a view row group by row group, keeping its schema and footer identity."""
    parquet = pq.ParquetFile(source)
    schema, compression = _stamped_schema(parquet.schema_arrow)
    codec = None if compression == "none" else compression
    with pq.ParquetWriter(str(destination), schema, compression=codec) as writer:
        for group in range(parquet.num_row_groups):
            table = parquet.read_row_group(group)
            writer.write_table(fill_batch(table).cast(schema))


def _fill_pg_table(table: pa.Table, fasta: FastaSequences, protein_peptides, report: ProteinPropertiesReport) -> pa.Table:
    names = table.schema.names
    report.pg_rows += table.num_rows
    anchors = table.column("anchor_protein").to_pylist()
    decoys = table.column("is_decoy").to_pylist() if "is_decoy" in names else [False] * table.num_rows
    coverage = table.column("sequence_coverage").to_pylist() if "sequence_coverage" in names else None
    weight = table.column("molecular_weight").to_pylist() if "molecular_weight" in names else None
    if coverage is None and weight is None:
        return table
    cache: dict[str, tuple[float | None, float | None]] = {}
    for row, anchor in enumerate(anchors):
        needs_coverage = coverage is not None and coverage[row] is None
        needs_weight = weight is not None and weight[row] is None
        if decoys[row] or not anchor or not (needs_coverage or needs_weight):
            continue
        report.pg_rows_eligible += 1
        if anchor not in cache:
            sequence = fasta.get(anchor)
            if sequence is None:
                cache[anchor] = (None, None)
            else:
                spans = [span for peptide in protein_peptides.get(anchor, ()) for span in peptide_occurrences(peptide, sequence)]
                cache[anchor] = (sequence_coverage_percent(sequence, spans), average_molecular_weight_kda(sequence))
        if fasta.get(anchor) is None:
            report.pg_anchors_not_in_fasta += 1
            _note_unmatched(report, anchor)
            continue
        cov, mw = cache[anchor]
        if needs_coverage and cov is not None:
            coverage[row] = cov
            report.pg_coverage_filled += 1
        if needs_weight and mw is not None:
            weight[row] = mw
            report.pg_molecular_weight_filled += 1
    if coverage is not None:
        index = names.index("sequence_coverage")
        table = table.set_column(index, table.schema.field(index), pa.array(coverage, type=table.schema.field(index).type))
    if weight is not None:
        index = names.index("molecular_weight")
        table = table.set_column(index, table.schema.field(index), pa.array(weight, type=table.schema.field(index).type))
    return table


def _fill_feature_table(table: pa.Table, fasta: FastaSequences, report: ProteinPropertiesReport) -> pa.Table:
    names = table.schema.names
    report.feature_rows += table.num_rows
    if "pg_positions" not in names or "pg_accessions" not in names:
        return table
    positions_column = table.column("pg_positions")
    if positions_column.null_count == 0:
        return table
    positions = positions_column.to_pylist()
    sequences = table.column("sequence").to_pylist()
    decoys = table.column("is_decoy").to_pylist() if "is_decoy" in names else [False] * table.num_rows
    groups = _accession_lists(table.column("pg_accessions"))
    cache: dict[tuple, list[dict] | None] = {}
    for row, existing in enumerate(positions):
        group = groups[row]
        if existing is not None or decoys[row] or not sequences[row] or not group:
            continue
        report.feature_rows_eligible += 1
        key = (sequences[row], group)
        if key not in cache:
            found: list[dict] = []
            for accession in group:
                protein = fasta.get(accession)
                for start, end in peptide_occurrences(sequences[row], protein):
                    found.append({"protein_accession": accession, "start": start, "end": end})
            cache[key] = found or None
        if cache[key] is None:
            report.feature_rows_without_match += 1
            continue
        positions[row] = cache[key]
        report.feature_positions_filled += 1
    index = names.index("pg_positions")
    return table.set_column(index, table.schema.field(index), pa.array(positions, type=table.schema.field(index).type))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _provenance_step(report: ProteinPropertiesReport, fasta_path: Path, views: list[str], step_order: int) -> dict:
    from qpx._version import __version__

    parameters = {
        "fasta": fasta_path.name,
        "fasta_sha256": _sha256(fasta_path),
        "fill_policy": "null values on target rows only; producer values are never overwritten",
        "pg_sequence_coverage_filled": report.pg_coverage_filled,
        "pg_molecular_weight_filled": report.pg_molecular_weight_filled,
        "pg_anchors_not_in_fasta": report.pg_anchors_not_in_fasta,
        "feature_pg_positions_filled": report.feature_positions_filled,
        "feature_rows_without_match": report.feature_rows_without_match,
    }
    return {
        "step_order": step_order,
        "step_category": "annotation",
        "step_name": "protein_properties_from_fasta",
        "tool_name": "qpx",
        "tool_version": __version__,
        "parameters": [{"key": key, "value": str(value)} for key, value in parameters.items()],
        "output_views": views,
    }


def _write_provenance(source: Path | None, destination: Path, step: dict) -> None:
    from qpx.writers.provenance import ProvenanceWriter

    steps = pq.read_table(source).to_pylist() if source is not None and source.is_file() else []
    step["step_order"] = max((s.get("step_order") or 0 for s in steps), default=0) + 1
    with ProvenanceWriter(destination, creator="qpx") as writer:
        writer.write_batch([*steps, step])


def annotate_protein_properties(
    dataset_dir: str | Path,
    prefix: str,
    fasta_path: str | Path,
    staging: str | Path,
) -> tuple[ProteinPropertiesReport, list[str]]:
    """Write FASTA-annotated pg/feature (and provenance) views for a dataset into ``staging``.

    Returns the report and the file names written, which the caller moves into
    place. Nothing in ``dataset_dir`` is modified here.
    """
    dataset_dir, staging, fasta_path = Path(dataset_dir), Path(staging), Path(fasta_path)
    fasta = FastaSequences.from_path(fasta_path)
    report = ProteinPropertiesReport(fasta_entries=fasta.entries, fasta_decoy_entries=fasta.decoy_entries)
    written: list[str] = []
    views: list[str] = []

    pg_path = dataset_dir / f"{prefix}.pg.parquet"
    feature_path = dataset_dir / f"{prefix}.feature.parquet"
    candidates = _peptide_protein_candidates(feature_path) if feature_path.is_file() else {}

    if pg_path.is_file():
        protein_peptides = _protein_peptides(candidates)
        name = pg_path.name
        _rewrite_view(pg_path, staging / name, lambda t: _fill_pg_table(t, fasta, protein_peptides, report))
        written.append(name)
        views.append("pg")
    if feature_path.is_file():
        name = feature_path.name
        _rewrite_view(feature_path, staging / name, lambda t: _fill_feature_table(t, fasta, report))
        written.append(name)
        views.append("feature")

    if written:
        provenance = dataset_dir / f"{prefix}.provenance.parquet"
        name = provenance.name
        _write_provenance(provenance, staging / name, _provenance_step(report, fasta_path, views, 0))
        written.append(name)

    match_rate = report.pg_match_rate
    if match_rate is not None and match_rate < 0.5:
        logger.warning(
            "Only %.1f%% of target protein groups were found in %s; is this the FASTA used for the search? "
            "Unmatched examples: %s",
            100 * match_rate,
            fasta_path.name,
            report.unmatched_examples,
        )
    return report, written
