"""Tests for mzML -> mz.parquet conversion (full spectra).

Covers ``SpectraMappingTransform.write_mz_parquet_from_dir``, the directory
entry point used by ``qpxc convert mz``. Mini mzML files are generated with
pyopenms so the tests are self-contained and fast.
"""

import gzip
import shutil

import pyarrow.parquet as pq
import pytest

pytest.importorskip("pyopenms")


def _make_mini_mzml(
    path, specs, native_id="controllerType=0 controllerNumber=1 scan={scan}", injection_time=None, ion_mobility=None
):
    """Write a minimal mzML at *path*.

    *specs* is a list of ``(ms_level, scan, rt_seconds)`` tuples. MS2+ spectra
    get a precursor. Native IDs default to the Thermo ``scan=N`` form so
    scan-number extraction can be exercised.
    """
    import pyopenms as oms

    exp = oms.MSExperiment()
    for ms_level, scan, rt in specs:
        spectrum = oms.MSSpectrum()
        spectrum.setMSLevel(ms_level)
        spectrum.setRT(rt)
        spectrum.setNativeID(native_id.format(scan=scan))
        if ms_level >= 2:
            precursor = oms.Precursor()
            precursor.setMZ(500.0 + scan)
            precursor.setCharge(2)
            spectrum.setPrecursors([precursor])
        spectrum.set_peaks(([100.0, 200.0, 300.0], [10.0, 20.0, 30.0]))
        exp.addSpectrum(spectrum)
    oms.MzMLFile().store(str(path), exp)
    if injection_time is not None:
        # Vendor converters write the ion injection time as a <scan> cvParam.
        cv_param = (
            f'<cvParam cvRef="MS" accession="MS:1000927" name="ion injection time" value="{injection_time}" '
            'unitAccession="UO:0000028" unitName="millisecond" unitCvRef="UO" />'
        )
        path.write_text(path.read_text().replace("<scan>", "<scan>" + cv_param))
    if ion_mobility is not None:
        # timsTOF converters (tdf2mzml) put the precursor's 1/K0 on the selected ion.
        cv_param = (
            f'<cvParam cvRef="MS" accession="MS:1002815" name="inverse reduced ion mobility" value="{ion_mobility}" '
            'unitCvRef="MS" unitAccession="MS:1002814" unitName="volt-second per square centimeter" />'
        )
        path.write_text(path.read_text().replace("<selectedIon>", "<selectedIon>" + cv_param))


@pytest.fixture
def mzml_dir(tmp_path):
    """Directory with one mini mzML: 1 MS1 + 2 MS2 spectra."""
    d = tmp_path / "mzml"
    d.mkdir()
    _make_mini_mzml(d / "run_a.mzML", [(1, 1, 60.0), (2, 2, 60.5), (2, 3, 61.0)])
    return d


def test_from_dir_all_levels(mzml_dir, tmp_path):
    from qpx.transforms.spectra_mapping import SpectraMappingTransform

    out = tmp_path / "out.mz.parquet"
    SpectraMappingTransform(mzml_directory=mzml_dir).write_mz_parquet_from_dir(out)

    table = pq.read_table(str(out))
    assert table.num_rows == 3  # 1 MS1 + 2 MS2
    assert {"run_file_name", "scan", "ms_level", "mz", "intensity"} <= set(table.column_names)
    assert set(table.column("run_file_name").to_pylist()) == {"run_a"}
    assert sorted(table.column("scan").to_pylist()) == [[1], [2], [3]]
    # every spectrum keeps its peaks
    assert all(len(mz) == 3 for mz in table.column("mz").to_pylist())


def test_from_dir_ms2_only(mzml_dir, tmp_path):
    from qpx.transforms.spectra_mapping import SpectraMappingTransform

    out = tmp_path / "out_ms2.mz.parquet"
    SpectraMappingTransform(mzml_directory=mzml_dir).write_mz_parquet_from_dir(out, ms_levels=[2])

    table = pq.read_table(str(out))
    assert table.num_rows == 2
    assert set(table.column("ms_level").to_pylist()) == {2}
    precursors = table.column("precursors").to_pylist()
    assert all(p and p[0]["selected_ion_mz"] > 0 for p in precursors)


def test_from_dir_reads_gzip(tmp_path):
    """pyopenms reads .mzML.gz directly; run name strips the .gz extension."""
    d = tmp_path / "mzml_gz"
    d.mkdir()
    plain = d / "run_b.mzML"
    _make_mini_mzml(plain, [(2, 10, 60.0)])
    with open(plain, "rb") as fin, gzip.open(str(d / "run_b.mzML.gz"), "wb") as fout:
        shutil.copyfileobj(fin, fout)
    plain.unlink()  # leave only the .gz

    from qpx.transforms.spectra_mapping import SpectraMappingTransform

    out = tmp_path / "out_gz.mz.parquet"
    SpectraMappingTransform(mzml_directory=d).write_mz_parquet_from_dir(out)

    table = pq.read_table(str(out))
    assert table.num_rows == 1
    assert table.column("run_file_name").to_pylist() == ["run_b"]
    assert table.column("scan").to_pylist() == [[10]]


def test_from_dir_empty_raises(tmp_path):
    from qpx.transforms.spectra_mapping import SpectraMappingTransform

    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(FileNotFoundError):
        SpectraMappingTransform(mzml_directory=empty).write_mz_parquet_from_dir(tmp_path / "x.mz.parquet")


@pytest.mark.parametrize(
    ("native_id", "expected"),
    [("index={scan}", [[4]]), ("frame=120 scan={scan}", [[120, 4]])],
)
def test_scan_keeps_native_id_components(tmp_path, native_id, expected):
    """mz.scan follows the PSM/feature convention, including tdf2mzml ``index=N`` IDs."""
    from qpx.transforms.spectra_mapping import SpectraMappingTransform

    d = tmp_path / "mzml"
    d.mkdir()
    _make_mini_mzml(d / "run_c.mzML", [(2, 4, 60.0)], native_id=native_id)
    out = tmp_path / "out.mz.parquet"
    SpectraMappingTransform(mzml_directory=d).write_mz_parquet_from_dir(out)

    assert pq.read_table(str(out)).column("scan").to_pylist() == expected


def test_ion_injection_time_comes_from_the_scan_acquisition(tmp_path):
    """The mzML scan's ion injection time (MS:1000927) is kept instead of 0."""
    from qpx.transforms.spectra_mapping import SpectraMappingTransform

    d = tmp_path / "mzml"
    d.mkdir()
    _make_mini_mzml(d / "run_d.mzML", [(1, 1, 60.0)], injection_time=55.0)
    _make_mini_mzml(d / "run_e.mzML", [(1, 1, 60.0)])
    out = tmp_path / "out.mz.parquet"
    SpectraMappingTransform(mzml_directory=d).write_mz_parquet_from_dir(out)

    rows = pq.read_table(str(out), columns=["run_file_name", "ion_injection_time"]).to_pylist()
    assert {row["run_file_name"]: row["ion_injection_time"] for row in rows} == {"run_d": 55.0, "run_e": 0.0}


def test_inverse_ion_mobility_comes_from_the_selected_ion(tmp_path):
    """PASEF MS2 spectra keep their 1/K0; spectra without one stay null."""
    from qpx.transforms.spectra_mapping import SpectraMappingTransform

    d = tmp_path / "mzml"
    d.mkdir()
    _make_mini_mzml(d / "run_g.mzML", [(2, 1, 60.0)], native_id="index={scan}", ion_mobility=0.88)
    _make_mini_mzml(d / "run_h.mzML", [(2, 1, 60.0)])
    out = tmp_path / "out.mz.parquet"
    SpectraMappingTransform(mzml_directory=d).write_mz_parquet_from_dir(out)

    rows = pq.read_table(str(out), columns=["run_file_name", "inverse_ion_mobility"]).to_pylist()
    mobility = {row["run_file_name"]: row["inverse_ion_mobility"] for row in rows}
    assert mobility["run_g"] == pytest.approx(0.88)
    assert mobility["run_h"] is None


def test_peak_arrays_use_byte_stream_split(tmp_path):
    """Spectrum m/z and intensity leaves are written with BYTE_STREAM_SPLIT."""
    from qpx.transforms.spectra_mapping import SpectraMappingTransform

    d = tmp_path / "mzml"
    d.mkdir()
    _make_mini_mzml(d / "run_f.mzML", [(1, 1, 60.0), (2, 2, 60.5)])
    out = tmp_path / "out.mz.parquet"
    SpectraMappingTransform(mzml_directory=d).write_mz_parquet_from_dir(out)

    row_group = pq.ParquetFile(str(out)).metadata.row_group(0)
    encodings = {row_group.column(i).path_in_schema: row_group.column(i).encodings for i in range(row_group.num_columns)}
    assert "BYTE_STREAM_SPLIT" in encodings["mz.list.element"]
    assert "BYTE_STREAM_SPLIT" in encodings["intensity.list.element"]
