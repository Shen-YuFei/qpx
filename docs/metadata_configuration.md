# Metadata Configuration Files

## Overview

The metadata configuration file feature generates comprehensive documentation of column mappings between source proteomics formats (MaxQuant, DIANN, quantms) and the QPX format. These CSV files serve as essential references for understanding data transformations and implementing software integrations.

### Purpose

The metadata configuration files document:

- **Column mappings**: Source format → QPX format transformations
- **Ontology terms**: Standard MS ontology associations (when available) 
- **Computation methods**: Data transformation and processing logic
- **Source information**: Original file context and workflow specifics

## File Format

The metadata configuration file is in **CSV format** with the following columns:

```csv
model_view,filename,column_name,ontology_accession,full_name,settings
```

### Fields

- **model_view**: Data model type (`psm`, `feature`, `pg`, `ibaq`)
- **filename**: Source filename where the data originates
- **column_name**: Column name in QPX format
- **ontology_accession**: Ontology term (e.g., `MS:1000041`, `UNIMOD:35`)
- **full_name**: Full descriptive name of the column
- **settings**: Description of how the value is computed or transformed

### Example

```csv
model_view,filename,column_name,ontology_accession,full_name,settings
psm,msms.txt,charge,MS:1000041,charge state,mapped from source column 'Charge'
psm,msms.txt,calculated_mz,,theoretical peptide mass-to-charge ratio,theoretical m/z calculated from sequence and modifications using PyOpenMS
feature,evidence.txt,intensities,,The intensity-based abundance of the peptide,structured array with sample_accession, channel, intensity
pg,proteinGroups.txt,additional_intensities,,Additional intensity values,includes LFQ intensity (issue #129: total_sum_unique_peptides, total_intensity_all_peptides, top3_intensity) and iBAQ values
```

## Supported Workflows

The metadata configuration generator supports five proteomics workflows, implementing the standardization requirements from [issue #129](https://github.com/bigbio/qpx/issues/129):

| Workflow | Description | Primary Use Case |
|----------|-------------|------------------|
| **maxquant** | MaxQuant proteomics analysis | Complete DDA proteomics pipeline |
| **diann** | DIA-NN Data-Independent Acquisition | DIA proteomics analysis |  
| **quantms-lfq** | Label-Free Quantification | Quantitative proteomics (LFQ) |
| **quantms-tmt** | Tandem Mass Tags quantification | Multiplexed quantitative proteomics |
| **quantms-psm** | Identification-focused pipeline | Protein identification studies |

## Usage

### Command Line Interface

Generate metadata configuration for a specific workflow:

```bash
# Generate for MaxQuant
qpxc utils metadata generate --workflow maxquant --output maxquant_config.csv

# Generate for DIANN
qpxc utils metadata generate --workflow diann --output diann_config.csv

# Generate for quantms-LFQ
qpxc utils metadata generate --workflow quantms-lfq --output quantms_lfq_config.csv

# Generate for quantms-TMT
qpxc utils metadata generate --workflow quantms-tmt --output quantms_tmt_config.csv

# Generate only for specific model (psm, feature, or pg)
qpxc utils metadata generate --workflow maxquant --model pg --output maxquant_pg_config.csv
```

### Python API

```python
from qpx.core.metadata import WorkflowMetadataGenerator

# Create generator for MaxQuant workflow
generator = WorkflowMetadataGenerator(workflow="maxquant")

# Generate metadata for all models (PSM, Feature, PG)
generator.generate_all_metadata()

# Save to CSV file
generator.generate_file("maxquant_metadata.csv")

# Or generate for specific model only
generator.clear()
generator.generate_pg_metadata()
generator.generate_file("maxquant_pg_only.csv")
```

### Custom Metadata Generation

For custom workflows or additional metadata:

```python
from qpx.core.metadata import MetadataSchemaGenerator

generator = MetadataSchemaGenerator()

# Add custom column metadata
generator.add_column_metadata(
    model_view="psm",
    filename="custom_file.txt",
    column_name="custom_score",
    ontology_accession="MS:9999999",
    full_name="Custom scoring metric",
    settings="computed using proprietary algorithm"
)

# Generate CSV file
generator.generate_file("custom_metadata.csv")
```

## Protein Quantification Standardization (Issue #129)

The metadata configuration files document the protein quantification methods defined in [issue #129](https://github.com/bigbio/qpx/issues/129):

### MaxQuant
- **LFQ intensity**: `total_sum_unique_peptides`, `total_intensity_all_peptides`, `top3_intensity`
- **iBAQ**: Intensity-based absolute quantification values
- **RAW intensity**: Raw protein intensity values

### DIANN
- **PG.MaxLFQ**: MaxLFQ algorithm output
- **PG.Quantity**: Direct quantity measurements
- Optional harmonized metrics: `total_all_peptides_intensity`, `top3_intensity`

### quantms-LFQ
- Default: `total_sum_unique_peptides`
- Additional: `total_intensity_all_peptides`, `top3_intensity`

### quantms-TMT
- Default: `top3_intensity` (from TMT reporter ions)
- Additional: `total_unique_peptides_intensity`, `total_all_peptides_intensity`

## Testing

Run the test suite to verify metadata generation:

```bash
# Run all metadata tests
python -m pytest tests/test_metadata_schema.py -v

# Run specific test
python -m pytest tests/test_metadata_schema.py::test_workflow_metadata_generator_maxquant -v
```

Test output files are saved to `configuration_test/` directory (automatically ignored by git).

## Example Script

See `examples/generate_metadata_config.py` for a complete example that generates metadata configuration files for all supported workflows.

```bash
# Run the example
python examples/generate_metadata_config.py
```

## Use Cases

### 1. Documentation
Use metadata configuration files to document your data processing pipeline and transformation logic.

### 2. Validation
Verify that column mappings and computations are correctly implemented by reviewing the generated metadata file.

### 3. Integration
Software tools that don't have access to OLS ontology terms can use the `settings` column to understand how data was computed.

### 4. Cross-workflow Comparison
Compare how different workflows handle protein quantification and other metrics by examining their respective metadata configuration files.

## Notes

- The `ontology_accession` field may be empty if no standard ontology term is available
- The `settings` column is particularly useful for documenting computed or transformed fields
- Metadata files are in standard CSV format (comma-separated)
- CSV files can be easily opened in Excel, Google Sheets, or any text editor
- All metadata is extracted from PyArrow schemas with additional workflow-specific mappings
