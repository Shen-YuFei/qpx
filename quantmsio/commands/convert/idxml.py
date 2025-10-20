import click
from pathlib import Path
from typing import Optional, Union
import glob
import logging
import tempfile
import shutil

from quantmsio.core.project import create_uuid_filename
from quantmsio.core.idxml import IdXML, merge_idxml_parquet_files


@click.command(
    "convert-idxml",
    short_help="Convert IdXML to PSM parquet file in quantms io",
)
@click.option(
    "--idxml-file",
    help="the IdXML file containing identifications",
    required=True,
)
@click.option(
    "--output-folder",
    help="Folder where the parquet file will be generated",
    required=True,
)
@click.option(
    "--mzml-file",
    help="Optional mzML to attach spectra by scan",
    required=False,
)
@click.option(
    "--output-prefix-file",
    help="Prefix of the parquet file needed to generate the file name",
    required=False,
)
@click.option(
    "--spectral-data",
    help="Spectral data fields (optional)",
    is_flag=True,
)
def convert_idxml_file(
    idxml_file: Union[Path, str],
    output_folder: str,
    mzml_file: Optional[Union[Path, str]],
    output_prefix_file: Optional[str],
    spectral_data: bool = False,
) -> None:

    if idxml_file is None or output_folder is None:
        raise click.UsageError("Please provide all the required parameters")

    if not output_prefix_file:
        output_prefix_file = "psm"

    parser = IdXML(
        idxml_path=idxml_file, mzml_path=mzml_file, spectral_data=spectral_data
    )
    output_path = (
        f"{output_folder}/{create_uuid_filename(output_prefix_file, '.psm.parquet')}"
    )
    parser.to_parquet(output_path)


@click.command(
    "convert-idxml-batch",
    short_help="Convert multiple IdXML files to a single merged PSM parquet file",
)
@click.option(
    "--idxml-folder",
    help="Folder containing IdXML files to convert",
    required=False,
    type=click.Path(exists=True, file_okay=False),
)
@click.option(
    "--idxml-files",
    help="Comma-separated list of IdXML file paths",
    required=False,
)
@click.option(
    "--output-folder",
    help="Folder where the merged parquet file will be generated",
    required=True,
)
@click.option(
    "--output-prefix-file",
    help="Prefix of the parquet file needed to generate the file name",
    required=False,
    default="merged-psm",
)
@click.option(
    "--mzml-folder",
    help="Optional folder containing mzML files to attach spectra by scan",
    required=False,
    type=click.Path(exists=True, file_okay=False),
)
@click.option(
    "--mzml-files",
    help="Comma-separated list of mzML file paths",
    required=False,
)
@click.option(
    "--verbose",
    help="Enable verbose logging",
    is_flag=True,
)
def convert_idxml_batch(
    idxml_folder: Optional[str],
    idxml_files: Optional[str],
    output_folder: str,
    output_prefix_file: str,
    mzml_folder: Optional[str],
    mzml_files: Optional[str],
    verbose: bool,
) -> None:
    """
    Convert multiple IdXML files to a single merged PSM parquet file.

    Input: --idxml-folder (all files) or --idxml-files (specific files)
    Optional: --mzml-folder or --mzml-files (for spectral data)
    Matching: Folder mode uses filename, file list mode uses position/index
    """
    # Setup logging
    logger = logging.getLogger(__name__)
    if verbose:
        logging.basicConfig(
            level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
        )
    else:
        logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")

    # Validate input parameters
    if not idxml_folder and not idxml_files:
        raise click.UsageError("Please provide either --idxml-folder or --idxml-files")

    if idxml_folder and idxml_files:
        raise click.UsageError(
            "Please provide only one of --idxml-folder or --idxml-files"
        )

    if mzml_folder and mzml_files:
        raise click.UsageError(
            "Please provide only one of --mzml-folder or --mzml-files"
        )

    # Get list of idXML files
    idxml_file_paths = []
    if idxml_folder:
        idxml_pattern = str(Path(idxml_folder) / "*.idXML")
        idxml_file_paths = glob.glob(idxml_pattern)
        if not idxml_file_paths:
            raise click.UsageError(f"No .idXML files found in folder: {idxml_folder}")
        logger.info(f"Found {len(idxml_file_paths)} idXML files in {idxml_folder}")
    else:
        idxml_file_paths = [f.strip() for f in idxml_files.split(",")]
        # Validate all files exist
        for idxml_file in idxml_file_paths:
            if not Path(idxml_file).exists():
                raise click.UsageError(f"IdXML file not found: {idxml_file}")
        logger.info(f"Processing {len(idxml_file_paths)} specified idXML files")

    # Get list of mzML files if provided
    mzml_file_paths = []
    mzml_map = None
    use_index_mapping = False  # Flag to indicate if we use index-based mapping

    if mzml_folder:
        # Will search for matching files later by basename
        pass
    elif mzml_files:
        mzml_file_paths = [f.strip() for f in mzml_files.split(",")]
        # Validate all files exist
        for mzml_file in mzml_file_paths:
            if not Path(mzml_file).exists():
                raise click.UsageError(f"mzML file not found: {mzml_file}")
        logger.info(f"Using {len(mzml_file_paths)} specified mzML files")

        # When using file lists, map by index (position) instead of basename
        if idxml_files:  # Both are file lists
            use_index_mapping = True
            logger.info(
                "Using index-based mapping: mzML files will be matched by position in the list"
            )
        else:  # idxml is folder, mzml is list - use basename matching
            mzml_map = {Path(f).stem: f for f in mzml_file_paths}
            logger.info("Using basename matching for mzML files")

    # Create output folder if it doesn't exist
    Path(output_folder).mkdir(parents=True, exist_ok=True)

    # Create a temporary directory for individual parquet files
    temp_dir = tempfile.mkdtemp(prefix="idxml_batch_")
    logger.info(f"Using temporary directory: {temp_dir}")

    try:
        # Convert each idXML file to parquet
        parquet_files = []
        for i, idxml_file in enumerate(idxml_file_paths):
            logger.info(f"Processing {i+1}/{len(idxml_file_paths)}: {idxml_file}")

            # Find corresponding mzML file
            mzml_file = None

            if use_index_mapping:
                # Index-based mapping: use position in list (no filename requirement)
                if i < len(mzml_file_paths):
                    mzml_file = mzml_file_paths[i]
                    logger.info(f"Matched by index [{i}]: {Path(mzml_file).name}")
                else:
                    logger.warning(f"No mzML file at index {i} for {idxml_file}")
            elif mzml_folder:
                # Basename matching in folder
                idxml_basename = Path(idxml_file).stem
                mzml_pattern = str(Path(mzml_folder) / f"{idxml_basename}.mzML")
                mzml_matches = glob.glob(mzml_pattern)
                if mzml_matches:
                    mzml_file = mzml_matches[0]
                    logger.info(f"Matched by basename: {Path(mzml_file).name}")
                else:
                    logger.warning(
                        f"No matching mzML file found for {idxml_file} (basename: {idxml_basename})"
                    )
            elif mzml_map:
                # Basename matching in file list (when idxml is folder, mzml is list)
                idxml_basename = Path(idxml_file).stem
                if idxml_basename in mzml_map:
                    mzml_file = mzml_map[idxml_basename]
                    logger.info(f"Matched by basename: {Path(mzml_file).name}")
                else:
                    logger.warning(
                        f"No matching mzML file found for {idxml_file} (basename: {idxml_basename})"
                    )

            # Parse and convert idXML to parquet
            try:
                # If mzML file is provided, automatically enable spectral data
                spectral_data = mzml_file is not None

                parser = IdXML(
                    idxml_path=idxml_file,
                    mzml_path=mzml_file,
                    spectral_data=spectral_data,
                )

                # Generate temporary parquet file name
                temp_parquet_file = (
                    Path(temp_dir)
                    / f"temp_{i}_{create_uuid_filename('psm', '.psm.parquet')}"
                )
                parser.to_parquet(str(temp_parquet_file))
                parquet_files.append(str(temp_parquet_file))

                logger.info(
                    f"Converted {idxml_file} -> {temp_parquet_file} ({parser.get_psm_count()} PSMs)"
                )
            except Exception as e:
                logger.error(f"Failed to convert {idxml_file}: {e}")
                if verbose:
                    raise
                continue

        # Merge all parquet files
        if not parquet_files:
            raise click.ClickException("No parquet files were successfully generated")

        logger.info(f"Merging {len(parquet_files)} parquet files...")
        output_path = Path(output_folder) / create_uuid_filename(
            output_prefix_file, ".psm.parquet"
        )
        merge_idxml_parquet_files(parquet_files, str(output_path))

        logger.info(f"Successfully created merged parquet file: {output_path}")

    finally:
        # Clean up temporary directory
        shutil.rmtree(temp_dir)
        logger.info(f"Cleaned up temporary directory: {temp_dir}")
