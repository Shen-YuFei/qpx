"""SQL safety tests for MuData label detection."""

from unittest.mock import Mock

import pytest

from qpx.mudata import _is_multiplexed


@pytest.mark.parametrize("table", ["feature", "pg"])
def test_is_multiplexed_uses_static_query_whitelist(table):
    """Supported tables select a fixed SQL statement."""
    engine = Mock()
    engine.execute.return_value.fetchall.return_value = [(2,)]

    assert _is_multiplexed(engine, table)
    engine.execute.assert_called_once()
    query = engine.execute.call_args.args[0]
    assert f"FROM {table}" in query
    if table == "pg":
        # pg is flattened since 1.1: a scalar ``label`` column (no intensities
        # list to UNNEST).
        assert "label" in query
    else:
        assert "i.label" in query


def test_is_multiplexed_rejects_unlisted_identifiers():
    """Unlisted identifiers are rejected before SQL execution."""
    engine = Mock()

    with pytest.raises(KeyError):
        _is_multiplexed(engine, "feature; DROP TABLE pg")
    engine.execute.assert_not_called()
