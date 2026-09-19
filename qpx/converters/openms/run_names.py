"""Normalize native OpenMS run references to the SDRF basename convention."""

from functools import cache

import pyarrow as pa

from qpx.core.constants import FEATURE, PG, PSM
from qpx.core.files import run_file_stem

_RUN_COLUMNS = {
    FEATURE: ("run_file_name", "id_run_file_name"),
    PSM: ("run_file_name",),
    PG: ("grouped_runs", "run_file_name"),
}


def normalize_run_names(table: pa.Table, view: str) -> pa.Table:
    """Normalize stored run references before deriving missing identities.

    Only acquisition-file suffixes are removed, so meaningful dots in run names
    survive. Missing references and grouped-run order are retained; this does
    not infer run membership or change supplied IDs and cross-references.
    """
    normalize = cache(run_file_stem)
    for name in _RUN_COLUMNS.get(view, ()):
        if name not in table.column_names:
            continue
        values = table.column(name).to_pylist()
        if name == "grouped_runs":
            values = [
                [normalize(run) if run is not None else None for run in group] if group is not None else None for group in values
            ]
        else:
            values = [normalize(run) if run is not None else None for run in values]
        field = table.schema.field(name)
        table = table.set_column(table.schema.get_field_index(name), field, pa.array(values, type=field.type))
    return table
