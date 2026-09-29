"""Regression coverage for source metadata preserved in QPX exports."""

import pyarrow.parquet as pq
import pytest

from qpx.converters.diann.converter import DiaNNConverter
from qpx.converters.sdrf import SdrfConverter


@pytest.mark.parametrize(
    ("log_text", "expected"),
    [
        ("DIA-NN 1.8.1\n", "1.8.1"),
        ("\n\nDIA-NN 2.5.1 Academia\n", "2.5.1"),
        ("Starting analysis\r\n\r\nDIA-NN 2.5.1 Academia\r\n", "2.5.1"),
        ("DIA-NN 2.5.1 Academia\nDIA-NN 1.8.1\n", "2.5.1"),
        ("Analysis complete\n", None),
        ("", None),
    ],
)
def test_diann_log_version_reaches_dataset_and_provenance(tmp_path, log_text, expected):
    """Blank lines or a preamble must not hide the first DIA-NN version."""
    log_path = tmp_path / "diann.log"
    log_path.write_text(log_text, encoding="utf-8")
    converter = DiaNNConverter(report_path="report.tsv", diann_log=str(log_path))

    converter.write_dataset(tmp_path, prefix="test")
    converter.write_provenance(tmp_path, prefix="test")

    dataset = pq.read_table(tmp_path / "test.dataset.parquet").to_pylist()
    provenance = pq.read_table(tmp_path / "test.provenance.parquet").to_pylist()
    assert dataset[0]["software_version"] == expected
    assert next(row for row in provenance if row["tool_name"] == "DIA-NN")["tool_version"] == expected


@pytest.mark.parametrize("provide_path", [False, True])
def test_diann_absent_log_leaves_version_unknown(tmp_path, provide_path):
    """An optional or unreadable log must not prevent metadata export."""
    log_path = str(tmp_path / "missing.log") if provide_path else None
    converter = DiaNNConverter(report_path="report.tsv", diann_log=log_path)
    converter.write_dataset(tmp_path, prefix="test")

    assert pq.read_table(tmp_path / "test.dataset.parquet").to_pylist()[0]["software_version"] is None


@pytest.mark.parametrize(
    ("compounds", "expected"),
    [
        (("E.coli 200 ng", "not applicable"), "E.coli 200 ng"),
        (("not applicable", "E.coli 200 ng"), "E.coli 200 ng"),
        (("E.coli 200 ng", "yeast 100 ng"), "E.coli 200 ng; yeast 100 ng"),
        (("E.coli 200 ng", "E.coli 200 ng"), "E.coli 200 ng"),
        (("E.coli 200 ng", ""), "E.coli 200 ng"),
        (("not applicable", "not applicable"), "not applicable"),
        (("", ""), None),
    ],
)
def test_repeated_sdrf_characteristics_preserve_values(tmp_path, compounds, expected):
    """Merge real compounds across duplicate columns without leaking samples."""
    sdrf_path = tmp_path / "repeated.sdrf.tsv"
    first_row = "\t".join(("S1", "Homo sapiens", "cell", *compounds)) + "\n"
    sdrf_path.write_text(
        "source name\tcharacteristics[organism]\tcharacteristics[organism part]\t"
        "characteristics[spiked compound]\tcharacteristics[spiked compound]\n"
        + first_row * 2
        + "S2\tHomo sapiens\tcell\tE.coli 500 ng\tnot applicable\n",
        encoding="utf-8",
    )
    sample_path = tmp_path / "sample.parquet"
    with SdrfConverter(duckdb_threads=1) as converter:
        converter.convert(str(sdrf_path), sample_output=str(sample_path))

    records = pq.read_table(sample_path).to_pylist()
    assert {row["sample_accession"]: row["spiked_compound"] for row in records} == {
        "S1": expected,
        "S2": "E.coli 500 ng",
    }
