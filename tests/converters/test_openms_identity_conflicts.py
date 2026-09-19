"""Ambiguous native identities must not replace any existing core output."""

from copy import deepcopy

import pytest

from qpx.converters.openms.converter import OpenMSConverter
from qpx.writers import FeatureWriter, PgWriter, PsmWriter
from tests.conftest import make_feature_record, make_pg_record, make_psm_record

_VIEW_INPUTS = {
    "feature": (FeatureWriter, make_feature_record),
    "psm": (PsmWriter, make_psm_record),
    "pg": (PgWriter, make_pg_record),
}


def _write_source_bundle(folder, conflict_view, conflict):
    """Put the conflicting records in separate row groups of a full bundle."""
    for view, (writer_class, factory) in _VIEW_INPUTS.items():
        first = factory()
        records = [first]
        if view == conflict_view:
            if conflict == "provided_id":
                first["feature_id"] = 101
            second = deepcopy(first)
            if conflict == "provided_id":
                second["rt"] += 1.0
            elif conflict == "different_match":
                second.update(rt=121.5, posterior_error_probability=0.2)
            elif conflict == "different_quantity":
                second["intensities"][0]["intensity"] += 1.0
            records.append(second)
        kwargs = {"override_provided_ids": False} if view == "feature" else {}
        with writer_class(folder / f"native.{view}.parquet", batch_size=1, **kwargs) as writer:
            for record in records:
                writer.write_batch([record])


@pytest.mark.parametrize(
    ("view", "conflict"),
    [
        ("feature", "identical"),
        ("feature", "provided_id"),
        ("psm", "different_match"),
        ("pg", "different_quantity"),
    ],
)
def test_native_identity_conflict_preserves_entire_existing_bundle(tmp_path, view, conflict):
    """Do not deduplicate records or invent IDs when source identity is ambiguous."""
    source = tmp_path / "source"
    source.mkdir()
    _write_source_bundle(source, view, conflict)
    original_inputs = {path: path.read_bytes() for path in source.iterdir()}
    output = tmp_path / "output"
    output.mkdir()
    original_outputs = {}
    for core_view in _VIEW_INPUTS:
        path = output / f"openms.{core_view}.parquet"
        original_outputs[path] = f"existing {core_view} output".encode()
        path.write_bytes(original_outputs[path])

    with pytest.raises(ValueError, match="Cannot safely identify native OpenMS .*duplicate.*openms-consensus"):
        OpenMSConverter(source).convert(output)

    assert {path: path.read_bytes() for path in output.iterdir()} == original_outputs
    assert {path: path.read_bytes() for path in source.iterdir()} == original_inputs
    assert not list(output.rglob("*.tmp"))
