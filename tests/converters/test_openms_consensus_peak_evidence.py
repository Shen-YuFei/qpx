"""Regression coverage for confidence scopes and duplicate measured-peak evidence."""

import pytest
from defusedxml import ElementTree

from qpx.converters.openms_consensus import streaming as stream_reader
from tests.converters.test_openms_consensus_merged_ids import (
    _HEADER,
    _assert_valid,
    _convert,
    _element,
    _maps,
    _percolator_run,
    _percolator_xml,
    _rows,
    _score,
    _xcorr_pid,
)


def _scoped_confidence_xml(second_value, same_source):
    """The same peptide in two runs, analysed either together or independently."""
    source = "PI_0" if same_source else "PI_1"
    second = _xcorr_pid("PeptideIdentification", 2, "PEPTIDEK", second_value, second_value)
    second = second.replace("PI_0", source).replace('name="map_index" value="0"', 'name="map_index" value="1"')
    header = _HEADER + _percolator_run("PI_0", "run_A.mzML,run_B.mzML" if same_source else "run_A.mzML", True)
    if not same_source:
        header += _percolator_run("PI_1", "run_B.mzML", True)
    elements = _element(0, [(0, 100, 1000)], [_xcorr_pid("PeptideIdentification", 1, "PEPTIDEK", 0.001, 0.01)])
    elements += _element(1, [(1, 100, 2000)], [second])
    return (
        header + _maps(["run_A.mzML", "run_B.mzML"]) + f"<consensusElementList>{elements}</consensusElementList></consensusXML>"
    )


@pytest.mark.parametrize("streaming", [False, True])
@pytest.mark.parametrize("same_source", [False, True])
@pytest.mark.parametrize("second_value", [0.2, 1.0])
def test_confidence_stays_with_its_fdr_source(tmp_path, streaming, same_source, second_value):
    """Only a shared analysis may propagate confidence across runs; placeholders stay null otherwise."""
    source = tmp_path / "sources.consensusXML"
    source.write_text(_scoped_confidence_xml(second_value, same_source))
    out, written = _convert(tmp_path, str(source), streaming)
    features = {rec["run_file_name"]: rec for rec in _rows(written, "feature")}
    expected_q = 0.001 if same_source else (second_value if second_value < 1 else None)
    expected_pep = 0.01 if same_source else expected_q
    assert features["run_A"]["peptide_qvalue"] == pytest.approx(0.001)
    assert features["run_B"]["peptide_qvalue"] == pytest.approx(expected_q)
    assert features["run_B"]["posterior_error_probability"] == pytest.approx(expected_pep)
    psm = next(rec for rec in _rows(written, "psm") if rec["run_file_name"] == "run_B")
    assert _score(psm, "peptide_qvalue") == pytest.approx(expected_q)
    _assert_valid(out, ("feature", "psm"))


@pytest.mark.parametrize("peptide_level", [False, True])
@pytest.mark.parametrize("unassigned_last", [False, True])
def test_streaming_uses_one_context_pass_and_one_output_pass(tmp_path, monkeypatch, peptide_level, unassigned_last):
    """Recovered IDs and peptide confidence share a pre-pass, even with trailing unassigned IDs."""
    source = tmp_path / "passes.consensusXML"
    root = ElementTree.fromstring(_percolator_xml(peptide_level))
    if unassigned_last:
        pid = root.find("UnassignedPeptideIdentification")
        root.remove(pid)
        root.append(pid)
    source.write_text(ElementTree.tostring(root, encoding="unicode"))
    parse = stream_reader.iterparse
    completed = []

    def count_complete_parses(*args, **kwargs):
        yield from parse(*args, **kwargs)
        completed.append(True)

    monkeypatch.setattr(stream_reader, "iterparse", count_complete_parses)
    _, written = _convert(tmp_path, str(source), True)
    assert len(completed) == 2
    feature = next(rec for rec in _rows(written, "feature") if rec["peptidoform"] == "PEPTIDEK")
    assert feature["posterior_error_probability"] == pytest.approx(0.01 if peptide_level else 0.004)
