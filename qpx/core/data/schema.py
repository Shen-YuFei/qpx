"""Schema base classes for QPX data structures.

ViewSchema instances are loaded from YAML specs at import time.
FieldDef is a simple container for a resolved Arrow field.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field

import pyarrow as pa
import pyarrow.compute as pc

from qpx.core.sql import sql_build

# ---------------------------------------------------------------------------
# Validation result types
# ---------------------------------------------------------------------------


@dataclass
class ValidationIssue:
    """A single validation issue found during schema validation."""

    structure: str
    check: str  # "missing_column" | "type_mismatch" | "null_values" | "duplicate_pk"
    severity: str  # "error" | "warning"
    column: str | None
    message: str


@dataclass
class ValidationResult:
    """Result of validating a structure against its schema."""

    structure: str
    issues: list[ValidationIssue] = field(default_factory=list)

    @property
    def is_valid(self) -> bool:
        """True if no errors (warnings are acceptable)."""
        return not any(i.severity == "error" for i in self.issues)

    @property
    def errors(self) -> list[ValidationIssue]:
        return [i for i in self.issues if i.severity == "error"]

    @property
    def warnings(self) -> list[ValidationIssue]:
        return [i for i in self.issues if i.severity == "warning"]

    @property
    def summary(self) -> str:
        n_err = len(self.errors)
        n_warn = len(self.warnings)
        status = "VALID" if self.is_valid else "INVALID"
        parts = []
        if n_err:
            parts.append(f"{n_err} error{'s' if n_err != 1 else ''}")
        if n_warn:
            parts.append(f"{n_warn} warning{'s' if n_warn != 1 else ''}")
        detail = f" ({', '.join(parts)})" if parts else ""
        return f"{self.structure}: {status}{detail}"


# ---------------------------------------------------------------------------
# Field and schema definitions
# ---------------------------------------------------------------------------


class FieldDef:
    """A resolved schema field definition."""

    __slots__ = ("name", "arrow_type", "nullable", "optional", "doc")

    def __init__(
        self,
        name: str,
        arrow_type: pa.DataType,
        *,
        nullable: bool = True,
        optional: bool = False,
        doc: str = "",
    ):
        self.name = name
        self.arrow_type = arrow_type
        self.nullable = nullable
        self.optional = optional
        self.doc = doc

    def to_arrow_field(self) -> pa.Field:
        metadata = {"doc": self.doc} if self.doc else None
        return pa.field(self.name, self.arrow_type, nullable=self.nullable, metadata=metadata)


def _types_compatible(expected: pa.DataType, actual: pa.DataType) -> bool:
    """Check type compatibility, tolerating nullability differences in nested structs.

    Parquet does not preserve inner-struct field nullability, so we compare
    types structurally while ignoring the ``nullable`` flag on struct children.
    """
    if expected == actual:
        return True

    # list<X> — compare value types recursively
    if isinstance(expected, pa.ListType) and isinstance(actual, pa.ListType):
        return _types_compatible(expected.value_type, actual.value_type)

    # map<K, V>
    if isinstance(expected, pa.MapType) and isinstance(actual, pa.MapType):
        return _types_compatible(expected.key_type, actual.key_type) and _types_compatible(expected.item_type, actual.item_type)

    # struct — compare field names, types (ignoring nullable)
    if isinstance(expected, pa.StructType) and isinstance(actual, pa.StructType):
        if expected.num_fields != actual.num_fields:
            return False
        for i in range(expected.num_fields):
            ef = expected.field(i)
            af = actual.field(i)
            if ef.name != af.name:
                return False
            if not _types_compatible(ef.type, af.type):
                return False
        return True

    return False


def _canonicalize_primary_key_column(column: pa.ChunkedArray) -> pa.Array | pa.ChunkedArray:
    """Normalize list-valued keys to sorted JSON for stable, collision-safe comparison."""
    if not pa.types.is_list(column.type):
        return column

    values = [
        json.dumps(sorted(value, key=lambda item: (item is None, str(item)))) if value is not None else None
        for value in column.to_pylist()
    ]
    return pa.array(values, type=pa.string())


_SCAN_FORMATS = ("scan", "index", "nativeId")


def _declared_scan_format(schema: pa.Schema, structure: str) -> str | None:
    """The scan_format a PSM view declares in its footer; None when there is nothing to check."""
    if structure != "psm" or "scan" not in schema.names:
        return None
    declared = (schema.metadata or {}).get(b"scan_format")
    return declared.decode(errors="replace") if declared else None


def _scan_format_bounds(scan_format: str) -> tuple[int, int]:
    return (2, 4) if scan_format == "nativeId" else (1, 1)


def _scan_format_violations(scans: pa.Array | pa.ChunkedArray, scan_format: str) -> int:
    """Non-empty scan arrays whose component count the declared format does not allow."""
    if scan_format not in _SCAN_FORMATS or not (pa.types.is_list(scans.type) or pa.types.is_large_list(scans.type)):
        return 0  # A non-list scan is reported by the type check.
    lengths = pc.call_function("list_value_length", [scans])
    minimum, maximum = _scan_format_bounds(scan_format)
    outside = pc.call_function(
        "or",
        [pc.call_function("less", [lengths, pa.scalar(minimum)]), pc.call_function("greater", [lengths, pa.scalar(maximum)])],
    )
    invalid = pc.call_function("and", [pc.call_function("greater", [lengths, pa.scalar(0)]), outside])
    return pc.call_function("sum", [invalid]).as_py() or 0


def _scan_format_issues(scan_format: str | None, count: int, structure: str, severity: str) -> list[ValidationIssue]:
    """Check a declared PSM format without guessing legacy or Feature encodings."""
    if scan_format is None:
        return []
    if scan_format not in _SCAN_FORMATS:
        return [
            _issue(structure, "scan_format", "warning", "scan", f"Unknown scan_format {scan_format!r}; cardinality not checked")
        ]
    if not count:
        return []
    minimum, maximum = _scan_format_bounds(scan_format)
    message = f"{count} non-empty PSM scan array(s) violate scan_format={scan_format!r}: expected {minimum}–{maximum} components"
    return [_issue(structure, "scan_format", severity, "scan", message)]


def _query_pg_referential_issues(
    con,
    structure: str,
    severity: str,
) -> list[ValidationIssue]:
    """Query protein-group referential invariants from ``pg_validate``.

    Both are intrinsic to the pg table — no run→sample resolution:

    * **duplicate_grouped_run** — ``grouped_runs`` is a set of distinct raw files;
      it must contain no duplicate element.
    * **run_double_count** (run-disjointness) — within one
      ``(pg_accessions, label)`` — the protein-group *membership* that keys the
      pg_id, not merely the leader — the ``grouped_runs`` sets across rows must be
      disjoint, so every raw file contributes to at most one row and no
      measurement is counted twice. Keying on the full membership (canonicalized
      order-independently) matches the pg_id identity, so two genuinely distinct
      groups that share a leading protein are not conflated into a false
      double-count. This is the join-free form of the double-count check.
    """
    issues: list[ValidationIssue] = []
    dup_lists = con.execute(
        "SELECT COUNT(*) FROM pg_validate WHERE grouped_runs IS NOT NULL "
        "AND length(grouped_runs) <> length(list_distinct(grouped_runs))"
    ).fetchall()[0][0]
    if dup_lists:
        issues.append(
            ValidationIssue(
                structure=structure,
                check="duplicate_grouped_run",
                severity=severity,
                column="grouped_runs",
                message=f"{dup_lists} pg row(s) have duplicate raw files within grouped_runs (must be a set of distinct files)",
            )
        )
    double = con.execute(
        """
        WITH exploded AS (
            -- Deduplicate each row's own grouped_runs first (a within-row
            -- duplicate is already reported by duplicate_grouped_run); this
            -- way run_double_count measures only cross-row disjointness.
            -- Key on the canonical (order-independent) group MEMBERSHIP, so
            -- distinct groups sharing a leader are not conflated.
            SELECT list_sort(list_distinct(pg_accessions)) AS members, label,
                   UNNEST(list_distinct(grouped_runs)) AS run
            FROM pg_validate
        ),
        repeated AS (
            SELECT COUNT(*) AS c
            FROM exploded
            GROUP BY members, label, run
            HAVING COUNT(*) > 1
        )
        SELECT COALESCE(SUM(c - 1), 0) FROM repeated
        """
    ).fetchall()[0][0]
    if double:
        issues.append(
            ValidationIssue(
                structure=structure,
                check="run_double_count",
                severity=severity,
                column=None,
                message=(
                    f"{double} run occurrence(s) repeat across pg rows sharing the same "
                    "(pg_accessions, label): grouped_runs sets must be disjoint or protein "
                    "intensity is double-counted"
                ),
            )
        )
    return issues


def _pg_referential_issues(
    table: pa.Table,
    structure: str,
    severity: str,
) -> list[ValidationIssue]:
    """Validate PG referential invariants for an in-memory Arrow table."""
    names = set(table.schema.names)
    if not {"pg_accessions", "grouped_runs", "label"} <= names or len(table) == 0:
        return []

    import duckdb

    con = duckdb.connect()
    try:
        con.register("pg_validate", table)
        return _query_pg_referential_issues(con, structure, severity)
    finally:
        con.close()


def pg_referential_issues_from_parquet(
    path: str,
    structure: str,
    severity: str,
) -> list[ValidationIssue]:
    """Validate PG referential invariants against a file-backed Parquet relation."""
    import duckdb

    con = duckdb.connect()
    try:
        relation = con.read_parquet(path)
        if not {"pg_accessions", "grouped_runs", "label"} <= set(relation.columns):
            return []
        relation.create_view("pg_validate")
        return _query_pg_referential_issues(con, structure, severity)
    finally:
        con.close()


def pg_referential_issues_from_query(con, query: str, structure: str, severity: str) -> list[ValidationIssue]:
    """Validate PG referential invariants over the rows of *query* on *con*."""
    con.execute(sql_build("CREATE OR REPLACE TEMP VIEW pg_validate AS $query", query=query))
    try:
        return _query_pg_referential_issues(con, structure, severity)
    finally:
        con.execute("DROP VIEW pg_validate")


# Unit-interval columns of the site and proteoform views: q-values,
# probabilities and evidence fractions.
_UNIT_INTERVAL_COLUMNS = {
    "site": ("site_qvalue", "modified_probability", "evidence_fraction"),
    "proteoform": ("proteoform_qvalue",),
}
# Tolerance for a residue's evidence fractions summing to 1.
_FRACTION_SUM_TOLERANCE = 1e-6
# Closed vocabularies of the site and proteoform views. ``prior_evidence`` and
# ``source`` are open: producers name their own evidence classes and sources.
_VOCABULARIES = {
    "proteoform": {"kind": ("canonical", "isoform", "variant", "chain"), "evidence": ("specific", "shared")},
}


def _out_of_range_count(values: pa.ChunkedArray | pa.Array) -> int:
    """Non-null values that are not finite numbers in [0, 1]; 0 for non-float input."""
    if not pa.types.is_floating(values.type):
        return 0
    bad = pc.or_(
        pc.or_(pc.less(values, 0.0), pc.greater(values, 1.0)),
        pc.invert(pc.is_finite(values)),
    )
    return pc.sum(pc.cast(pc.fill_null(bad, False), pa.int64())).as_py() or 0


def _issue(structure: str, check: str, severity: str, column: str, message: str) -> ValidationIssue:
    return ValidationIssue(structure=structure, check=check, severity=severity, column=column, message=message)


def _unit_interval_issues(counts: dict[str, int], structure: str, severity: str) -> list[ValidationIssue]:
    """Flag q-values, probabilities and evidence fractions outside [0, 1] or not finite."""
    issues = []
    for name, bad in counts.items():
        if bad:
            message = f"{bad} value(s) of {name} are not finite numbers in [0, 1]"
            issues.append(_issue(structure, "value_out_of_range", severity, name, message))
    return issues


def _alternative_fractions(table: pa.Table) -> list[list | None] | None:
    """Each site row's alternative evidence fractions, or None when the column is absent or malformed."""
    if "alternatives" not in table.schema.names:
        return None
    alternatives = table.column("alternatives")
    if not (pa.types.is_list(alternatives.type) and pa.types.is_struct(alternatives.type.value_type)):
        return None
    index = alternatives.type.value_type.get_field_index("evidence_fraction")
    if index < 0 or not pa.types.is_floating(alternatives.type.value_type.field(index).type):
        return None
    return [None if alts is None else [a.get("evidence_fraction") for a in alts] for alts in alternatives.to_pylist()]


def _evidence_fraction_counts(table: pa.Table | pa.RecordBatch) -> tuple[int, int, int]:
    """Alternatives out of [0, 1], rows not summing to 1, and rows listed out of order.

    A residue's evidence is one pool: a row's fraction and its alternatives' sum to 1.
    Without alternatives, the row holds the whole pool, so its fraction must be 1
    (or null). Alternatives are listed largest first and are each in [0, 1].
    """
    fractions = _alternative_fractions(table)
    if fractions is None or "evidence_fraction" not in table.schema.names:
        return 0, 0, 0
    if not pa.types.is_floating(table.column("evidence_fraction").type):
        return 0, 0, 0
    own = table.column("evidence_fraction").to_pylist()
    out_of_range = unsummed = unordered = 0
    for row_own, alts in zip(own, fractions):
        alts = [a for a in (alts or []) if a is not None]
        out_of_range += sum(not (0.0 <= a <= 1.0) for a in alts)
        unordered += any(x < y for x, y in zip(alts, alts[1:]))
        if row_own is not None and abs(row_own + sum(alts) - 1.0) > _FRACTION_SUM_TOLERANCE:
            unsummed += 1
    return out_of_range, unsummed, unordered


def _evidence_fraction_issues(counts: tuple[int, int, int], structure: str, severity: str) -> list[ValidationIssue]:
    out_of_range, unsummed, unordered = counts
    issues = []
    if out_of_range:
        message = f"{out_of_range} alternative evidence fraction(s) are not finite numbers in [0, 1]"
        issues.append(_issue(structure, "value_out_of_range", severity, "alternatives", message))
    if unsummed:
        message = (
            f"{unsummed} row(s) have an evidence fraction and alternatives that do not sum to 1 "
            "(a row without alternatives holds the whole residue's evidence: fraction 1)"
        )
        issues.append(_issue(structure, "evidence_fraction_sum", severity, "evidence_fraction", message))
    if unordered:
        message = f"{unordered} row(s) list alternatives that are not in decreasing evidence fraction"
        issues.append(_issue(structure, "alternatives_order", severity, "alternatives", message))
    return issues


def _vocabulary_count(values: pa.Array | pa.ChunkedArray, allowed: tuple[str, ...]) -> int:
    """Non-null values outside a closed vocabulary; 0 for non-string input."""
    if not (pa.types.is_string(values.type) or pa.types.is_large_string(values.type)):
        return 0
    outside = pc.invert(pc.is_in(values, value_set=pa.array(allowed)))
    return pc.sum(pc.cast(pc.and_(outside, pc.is_valid(values)), pa.int64())).as_py() or 0


def _vocabulary_issues(counts: dict[str, int], structure: str, severity: str) -> list[ValidationIssue]:
    """Flag values outside a closed vocabulary (null is allowed)."""
    issues = []
    for name, bad in counts.items():
        if bad:
            message = f"{bad} value(s) of {name} are not one of {', '.join(_VOCABULARIES[structure][name])}"
            issues.append(_issue(structure, "value_not_in_vocabulary", severity, name, message))
    return issues


# Columns holding sequence edits, per view, and the keys one sequence is
# identified by within a row (pepmap: one per protein mapping).
_SEQUENCE_EDIT_COLUMNS = {"pepmap": ("sequence_edits", ("protein_accession",)), "proteoform": ("edits", ())}


def _is_float(t: pa.DataType) -> bool:
    return pa.types.is_floating(t)


def _is_int(t: pa.DataType) -> bool:
    return pa.types.is_integer(t)


def _is_str(t: pa.DataType) -> bool:
    return pa.types.is_string(t) or pa.types.is_large_string(t)


def _flattened_struct_column(table: pa.Table | pa.RecordBatch, column: str, fields: dict) -> pa.Table | None:
    """One row per list element of a list<struct> column, with its parent row index.

    ``fields`` maps each needed struct field to a type predicate. None when the
    column is absent or empty, not a list of structs, or a field is missing or of
    another type: the type check reports that, and semantic checks must not run
    arithmetic or string kernels on it.
    """
    if column not in table.schema.names or len(table) == 0:
        return None
    col = table.column(column)
    if not (pa.types.is_list(col.type) and pa.types.is_struct(col.type.value_type)):
        return None
    struct = col.type.value_type
    for name, ok in fields.items():
        index = struct.get_field_index(name)
        if index < 0 or not ok(struct.field(index).type):
            return None
    flat = pc.list_flatten(col)
    return pa.table({"row": pc.list_parent_indices(col), **{f: pc.struct_field(flat, f) for f in fields}})


def _modification_candidate_counts(table: pa.Table | pa.RecordBatch) -> tuple[int, int]:
    """Candidate probabilities out of [0, 1], and alternative groups whose probabilities sum above 1.

    Candidates without a group may co-occur (a label and an acetyl on one
    residue), so no sum is implied across them.
    """
    flat = _flattened_struct_column(
        table,
        "modification_candidates",
        {"position": _is_int, "probability": _is_float, "alternative_group": _is_int},
    )
    if flat is None or flat.num_rows == 0:
        return 0, 0
    bad = _out_of_range_count(flat.column("probability"))
    over = 0
    grouped = flat.filter(pc.is_valid(flat.column("alternative_group")))
    if grouped.num_rows:
        sums = grouped.group_by(["row", "alternative_group"]).aggregate([("probability", "sum")])
        over = pc.greater(sums.column("probability_sum"), 1.0 + _FRACTION_SUM_TOLERANCE)
        over = pc.sum(pc.cast(pc.fill_null(over, False), pa.int64())).as_py() or 0
    return bad, over


def _modification_candidate_issues(counts: tuple[int, int], structure: str, severity: str) -> list[ValidationIssue]:
    """Candidate probabilities are in [0, 1]; within an alternative group they sum to at most 1."""
    bad, over = counts
    issues = []
    if bad:
        message = f"{bad} modification candidate probabilities are not finite numbers in [0, 1]"
        issues.append(_issue(structure, "value_out_of_range", severity, "modification_candidates", message))
    if over:
        message = f"{over} alternative group(s) have candidate probabilities that sum to more than 1"
        issues.append(_issue(structure, "candidate_probability_sum", severity, "modification_candidates", message))
    return issues


def _malformed_edit_count(flat: pa.Table) -> int:
    """Edits that do not name a real reference span and residues, or change nothing."""
    start, end = flat.column("start"), flat.column("end")
    ref, alt = flat.column("reference"), flat.column("alternative")
    span = pc.add(pc.subtract(end, start), 1)
    letters = r"^[A-Z]*$"
    bad = pc.or_kleene(
        pc.or_kleene(pc.less(start, 1), pc.not_equal(pc.utf8_length(ref), span)),
        pc.or_kleene(
            pc.or_kleene(
                pc.invert(pc.match_substring_regex(ref, letters)),
                pc.invert(pc.match_substring_regex(alt, letters)),
            ),
            pc.equal(ref, alt),
        ),
    )
    return pc.sum(pc.cast(pc.fill_null(bad, True), pa.int64())).as_py() or 0


def _conflicting_edit_count(flat: pa.Table, keys: tuple[str, ...]) -> int:
    """Sequences whose edits are out of order or overlap, so they define no single sequence.

    An edit covers the reference interval ``[start, end]``; an insertion before
    residue ``s`` (``end = s - 1``) is the point ``s - 1/2``. Within one sequence
    (a row, or a row and protein mapping), each edit must start strictly after
    the previous one ends: two edits of one residue, overlapping spans and two
    insertions at one point all conflict.
    """
    previous_end: dict[tuple, float] = {}
    conflicting: set[tuple] = set()
    for record in flat.select(["row", *keys, "start", "end"]).to_pylist():
        start, end = record["start"], record["end"]
        if start is None or end is None:
            continue
        key = tuple(record[k] for k in ("row", *keys))
        low, high = (start - 0.5, start - 0.5) if end == start - 1 else (start, end)
        if key in previous_end and low <= previous_end[key]:
            conflicting.add(key)
        previous_end[key] = max(high, previous_end.get(key, high))
    return len(conflicting)


def _sequence_edit_counts(table: pa.Table | pa.RecordBatch, structure: str) -> tuple[int, int]:
    """Malformed sequence edits, and sequences whose edits overlap or are out of order.

    ``reference`` is the span ``start..end`` (1-based, inclusive) of the reference
    sequence, so its length is ``end - start + 1``; an insertion has an empty
    reference and ``end = start - 1``. Residues are one-letter codes, and an edit
    changes something.
    """
    column, keys = _SEQUENCE_EDIT_COLUMNS.get(structure, (None, ()))
    if column is None:
        return 0, 0
    fields = {"start": _is_int, "end": _is_int, "reference": _is_str, "alternative": _is_str}
    fields.update({k: _is_str for k in keys})
    flat = _flattened_struct_column(table, column, fields)
    if flat is None or flat.num_rows == 0:
        return 0, 0
    return _malformed_edit_count(flat), _conflicting_edit_count(flat, keys)


def _sequence_edit_issues(counts: tuple[int, int], structure: str, severity: str) -> list[ValidationIssue]:
    """Each sequence edit names a real reference span and residues, and the edits of
    one sequence are ordered and do not overlap.
    """
    malformed, conflicting = counts
    column = _SEQUENCE_EDIT_COLUMNS.get(structure, (None, ()))[0]
    issues = []
    if malformed:
        message = (
            f"{malformed} sequence edit(s) are malformed: start must be >= 1, reference must be the "
            "end - start + 1 residues of the span (empty with end = start - 1 for an insertion), "
            "residues must be one-letter codes, and an edit must change something"
        )
        issues.append(_issue(structure, "malformed_sequence_edit", severity, column, message))
    if conflicting:
        message = (
            f"{conflicting} sequence(s) have edits that overlap or are out of order, so they do not define "
            "one sequence: each edit must start after the previous one ends"
        )
        issues.append(_issue(structure, "conflicting_sequence_edits", severity, column, message))
    return issues


def _anchor_without_membership(table: pa.Table | pa.RecordBatch) -> int:
    """Count rows that carry an ``anchor_protein`` but no ``pg_accessions`` membership.

    The protein group (``pg_accessions``) is the analyte and the feature->pg join
    key (see docs/spec/pg.md -> Protein group semantics); ``anchor_protein`` is only
    a descriptive representative. A row that has been protein-mapped (non-null
    ``anchor_protein``) MUST therefore carry its full group membership, otherwise it
    cannot be attributed to a group. Applies to any view exposing both columns
    (``feature``, where both are nullable; ``pg`` already requires both).
    """
    names = table.schema.names
    if "anchor_protein" not in names or "pg_accessions" not in names or len(table) == 0:
        return 0
    anchor = table.column("anchor_protein")
    pg = table.column("pg_accessions")
    # A non-list pg_accessions is a type mismatch already reported by the type check
    # above; skip here rather than let list_value_length raise on malformed input.
    if not (pa.types.is_list(pg.type) or pa.types.is_large_list(pg.type)):
        return 0
    # Non-string anchor_protein is a type mismatch reported elsewhere.
    if not (pa.types.is_string(anchor.type) or pa.types.is_large_string(anchor.type)):
        return 0
    # list_value_length is null where pg_accessions is null; fill_null(0) folds a
    # null membership into "empty" so the boolean logic never propagates nulls.
    pg_len = pc.fill_null(pc.list_value_length(pg), 0)
    pg_missing = pc.equal(pg_len, 0)
    # Blank/whitespace anchors are unset (not a membership violation).
    anchor_set = pc.greater(pc.utf8_length(pc.utf8_trim_whitespace(pc.fill_null(anchor, ""))), 0)
    violation = pc.and_(anchor_set, pg_missing)
    return pc.sum(pc.cast(violation, pa.int64())).as_py() or 0


def _anchor_membership_issues(bad: int, structure: str, severity: str) -> list[ValidationIssue]:
    if not bad:
        return []
    message = (
        f"{bad} row(s) have anchor_protein set but no pg_accessions "
        "(protein-group membership). The group is the analyte and the "
        "feature->pg join key, so a row with a leading protein must carry "
        "its full pg_accessions membership."
    )
    return [_issue(structure, "anchor_without_membership", severity, "pg_accessions", message)]


def _primary_key_issues(
    keys: pa.Table | None,
    primary_key: tuple[str, ...],
    structure: str,
    severity: str,
) -> list[ValidationIssue]:
    """Return a duplicate-key issue for the primary-key columns of a whole view.

    *keys* is None when a primary-key column is missing.
    """
    if keys is None or keys.num_rows == 0:
        return []

    key_table = pa.table({column: _canonicalize_primary_key_column(keys.column(column)) for column in primary_key})
    total_count = len(key_table)
    unique_count = key_table.group_by(list(primary_key)).aggregate([]).num_rows
    if unique_count == total_count:
        return []

    duplicate_count = total_count - unique_count
    return [
        ValidationIssue(
            structure=structure,
            check="duplicate_pk",
            severity=severity,
            column=None,
            message=(
                f"Primary key ({', '.join(primary_key)}) has "
                f"{duplicate_count} duplicate row(s) "
                f"({unique_count} unique out of {total_count})"
            ),
        )
    ]


class _RowCounts:
    """Counts of the checks that read one row at a time, summed over a view's batches."""

    def __init__(self, structure: str, schema: pa.Schema, required: list[str]):
        names = schema.names
        self.structure = structure
        self.nulls = dict.fromkeys(required, 0)
        self.anchor = 0
        self.out_of_range = dict.fromkeys((name for name in _UNIT_INTERVAL_COLUMNS.get(structure, ()) if name in names), 0)
        self.fractions = (0, 0, 0)
        self.candidates = (0, 0)
        self.edits = (0, 0)
        self.vocabulary = dict.fromkeys((name for name in _VOCABULARIES.get(structure, {}) if name in names), 0)
        self.scan = 0

    def add(self, batch: pa.Table | pa.RecordBatch, scan_format: str | None) -> None:
        """Add one batch's counts."""
        for name in self.nulls:
            self.nulls[name] += batch.column(name).null_count
        self.anchor += _anchor_without_membership(batch)
        for name in self.out_of_range:
            self.out_of_range[name] += _out_of_range_count(batch.column(name))
        if self.structure == "site":
            self.fractions = tuple(a + b for a, b in zip(self.fractions, _evidence_fraction_counts(batch)))
        self.candidates = tuple(a + b for a, b in zip(self.candidates, _modification_candidate_counts(batch)))
        self.edits = tuple(a + b for a, b in zip(self.edits, _sequence_edit_counts(batch, self.structure)))
        for name in self.vocabulary:
            self.vocabulary[name] += _vocabulary_count(batch.column(name), _VOCABULARIES[self.structure][name])
        if scan_format is not None:
            self.scan += _scan_format_violations(batch.column("scan"), scan_format)

    def null_issues(self, rows: int, severity: str) -> list[ValidationIssue]:
        return [
            _issue(
                self.structure,
                "null_values",
                severity,
                name,
                f"Column '{name}' is non-nullable but contains {count} null value(s) out of {rows} rows",
            )
            for name, count in self.nulls.items()
            if count
        ]

    def value_issues(self, severity: str) -> list[ValidationIssue]:
        return [
            *_anchor_membership_issues(self.anchor, self.structure, severity),
            *_unit_interval_issues(self.out_of_range, self.structure, severity),
            *_evidence_fraction_issues(self.fractions, self.structure, severity),
            *_vocabulary_issues(self.vocabulary, self.structure, severity),
            *_modification_candidate_issues(self.candidates, self.structure, severity),
            *_sequence_edit_issues(self.edits, self.structure, severity),
        ]


class ViewSchema:
    """Schema for a QPX view, loaded from YAML."""

    def __init__(
        self,
        view_name: str,
        file_type: str,
        primary_key: list[str],
        fields: dict[str, FieldDef],
        extra_columns: bool = False,
        identity_composite: list[str] | None = None,
    ):
        self._view_name = view_name
        self._file_type = file_type
        self._primary_key = tuple(primary_key)
        self._fields = fields
        self._extra_columns = extra_columns
        self._identity_composite = tuple(identity_composite) if identity_composite else None
        self._arrow_schema_cache: pa.Schema | None = None

    @property
    def primary_key(self) -> tuple[str, ...]:
        """The view's primary-key column(s)."""
        return self._primary_key

    @property
    def view_name(self) -> str:
        """The schema view name."""
        return self._view_name

    @property
    def identity_composite(self) -> tuple[str, ...] | None:
        """Composite of existing fields the single-column identity id is derived
        from, as declared by the schema's ``identity_composite`` key. ``None``
        for views without a derived identity (mz, run, sample, ...)."""
        return self._identity_composite

    @property
    def __name__(self) -> str:
        """For compatibility with parametrized tests and repr."""
        parts = self._view_name.split("_")
        return "".join(p.capitalize() for p in parts) + "Schema"

    def __repr__(self) -> str:
        return f"ViewSchema({self._view_name!r}, fields={len(self._fields)})"

    def get_arrow_schema(self) -> pa.Schema:
        """Generate PyArrow schema from field definitions."""
        if self._arrow_schema_cache is None:
            arrow_fields = [fdef.to_arrow_field() for fdef in self._fields.values()]
            self._arrow_schema_cache = pa.schema(arrow_fields)
        return self._arrow_schema_cache

    def get_extended_arrow_schema(self, extra_column_names: list[str]) -> pa.Schema:
        """Return the base schema extended with additional string columns.

        Only allowed when ``extra_columns=True`` in the YAML spec.
        Extra columns are appended as nullable ``pa.string()`` fields.
        """
        if not self._extra_columns:
            raise ValueError(f"Schema '{self._view_name}' does not allow extra columns")
        base = self.get_arrow_schema()
        known = set(base.names)
        extra_fields = [pa.field(name, pa.string(), nullable=True) for name in extra_column_names if name not in known]
        if not extra_fields:
            return base
        return pa.schema(list(base) + extra_fields)

    def _iter_fields(self):
        """Yield (name, FieldDef) pairs in definition order."""
        yield from self._fields.items()

    # --- Validation --------------------------------------------------------

    def validate(self, table: pa.Table, *, strict: bool = False) -> list[str]:
        """Validate an Arrow table against this schema. Returns list of error strings.

        For structured results, use validate_full() instead.
        """
        result = self.validate_full(table, strict=strict)
        return [issue.message for issue in result.errors]

    def validate_full(self, table: pa.Table, *, strict: bool = False) -> ValidationResult:
        """Full schema validation returning structured results.

        Checks:
            1. Required columns are present
            2. Column types match the schema
            3. Non-nullable columns contain no null values
            4. Primary key is unique
            5. Non-empty PSM scans match a declared scan_format

        Args:
            table: The Arrow table to validate.
            strict: When True (used by ``qpxc validate`` and CI), nulls in
                non-nullable columns and duplicate primary keys are reported as
                ``error`` (so ``is_valid`` is False and the CLI exits non-zero).
                Defaults to False so that *writers* and *converters*, which
                persist source data as-produced, are not blocked by legitimate
                duplicate keys (e.g. a DIA-NN feature whose stripped-sequence PK
                repeats across modiforms) — strictness is an audit concern, not
                a write-time one. Missing columns and type mismatches are always
                errors regardless of ``strict``.
        """
        severity = "error" if strict else "warning"
        keys = table.select(list(self._primary_key)) if set(self._primary_key) <= set(table.schema.names) else None
        referential = _pg_referential_issues(table, self._view_name, severity) if self._view_name == "pg" else []
        return self.validate_batches(table.schema, [table], keys, referential, strict=strict)

    def validate_batches(
        self,
        schema: pa.Schema,
        batches,
        keys: pa.Table | None,
        referential: list[ValidationIssue],
        *,
        strict: bool = False,
    ) -> ValidationResult:
        """Validate a view read as Arrow batches sharing *schema*, with the checks of validate_full().

        The checks that read one row at a time are summed batch by batch, so the
        whole view is never held at once. The two that span rows are passed in for
        the whole view: *keys* holds its primary-key columns (None when one is
        missing) and *referential* its pg referential issues.
        """
        result = ValidationResult(structure=self._view_name)
        severity = "error" if strict else "warning"
        result.issues.extend(self._column_issues(schema))

        required = [name for name, fdef in self._fields.items() if name in schema.names and not fdef.nullable]
        counts = _RowCounts(self._view_name, schema, required)
        scan_format = _declared_scan_format(schema, self._view_name)
        rows = 0
        for batch in batches:
            rows += batch.num_rows
            counts.add(batch, scan_format)

        result.issues.extend(counts.null_issues(rows, severity))
        result.issues.extend(_primary_key_issues(keys, self._primary_key, self._view_name, severity))
        result.issues.extend(referential)
        result.issues.extend(counts.value_issues(severity))
        result.issues.extend(_scan_format_issues(scan_format, counts.scan, self._view_name, severity))
        return result

    def _column_issues(self, schema: pa.Schema) -> list[ValidationIssue]:
        """Required columns that are missing and columns whose type does not match."""
        issues = []
        for expected_field in self.get_arrow_schema():
            name = expected_field.name
            if name not in schema.names:
                if not self._is_optional(name):
                    issues.append(_issue(self._view_name, "missing_column", "error", name, f"Missing required column: {name}"))
                continue
            actual = schema.field(name)
            if not _types_compatible(expected_field.type, actual.type):
                message = f"Column '{name}': expected {expected_field.type}, got {actual.type}"
                issues.append(_issue(self._view_name, "type_mismatch", "error", name, message))
        return issues

    def _is_optional(self, field_name: str) -> bool:
        fdef = self._fields.get(field_name)
        return fdef is not None and fdef.optional
