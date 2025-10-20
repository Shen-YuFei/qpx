# Visualization Commands

Create various data visualization plots from quantms.io data.

## Overview

The `visualize` command group provides tools for creating publication-quality visualizations from quantms.io parquet files. All plots are saved in vector formats (SVG, PDF) for high-resolution output.

## Available Commands

All visualization commands are accessed through the `plot` subcommand:

- [psm-peptides](#psm-peptides) - Plot peptides by condition in LFQ
- [ibaq-distribution](#ibaq-distribution) - Plot iBAQ distribution
- [kde-intensity](#kde-intensity) - Plot KDE intensity distribution
- [peptide-distribution](#peptide-distribution) - Plot peptide distribution across proteins
- [box-intensity](#box-intensity) - Plot intensity box plots

---

## psm-peptides

Plot peptides by condition in label-free quantification (LFQ) experiments.

### Description

Creates a visualization showing the distribution of identified peptides across different experimental conditions. This plot helps assess data quality and completeness across samples.

### Parameters

| Parameter            | Type | Required | Default | Description                        |
| -------------------- | ---- | -------- | ------- | ---------------------------------- |
| `--psm-parquet-path` | Path | Yes      | -       | PSM parquet file path              |
| `--sdrf-path`        | Path | Yes      | -       | SDRF file path for metadata        |
| `--save-path`        | Path | Yes      | -       | Output image path (e.g., plot.svg) |

### Usage Examples

#### Basic Example

```bash
quantmsioc visualize plot psm-peptides \
    --psm-parquet-path tests/examples/parquet/psm.parquet \
    --sdrf-path tests/examples/quantms/dda-lfq-full/PXD007683-LFQ.sdrf.tsv \
    --save-path ./plots/peptides_by_condition.svg
```

#### Generate PDF Output

```bash
quantmsioc visualize plot psm-peptides \
    --psm-parquet-path ./output/psm.parquet \
    --sdrf-path ./metadata.sdrf.tsv \
    --save-path ./plots/peptides_by_condition.pdf
```

### Output

- **Format**: SVG or PDF (based on file extension in `--save-path`)
- **Content**: Bar plot showing peptide counts per condition
- **Axes**:
  - X-axis: Experimental conditions
  - Y-axis: Number of identified peptides

### Interpretation

- **High variation**: May indicate batch effects or quality issues
- **Low counts**: May suggest technical problems with specific samples
- **Consistent counts**: Indicates good data quality and reproducibility

### Best Practices

- Use SVG format for publications (scalable vector graphics)
- Verify SDRF metadata correctly defines experimental conditions
- Compare with expected peptide yields for your sample type

---

## ibaq-distribution

Plot the distribution of iBAQ (intensity-Based Absolute Quantification) values.

### Description

Creates a histogram or density plot showing the distribution of iBAQ values across proteins. Useful for quality control and understanding the dynamic range of protein quantification.

### Parameters

| Parameter         | Type   | Required | Default | Description                          |
| ----------------- | ------ | -------- | ------- | ------------------------------------ |
| `--ibaq-path`     | Path   | Yes      | -       | iBAQ file path                       |
| `--save-path`     | Path   | Yes      | -       | Output image path (e.g., plot.svg)   |
| `--select-column` | String | No       | -       | Specific column in iBAQ file to plot |

### Usage Examples

#### Plot All Samples

```bash
quantmsioc visualize plot ibaq-distribution \
    --ibaq-path tests/examples/AE/PXD016999.1-ibaq.tsv \
    --save-path ./plots/ibaq_distribution.svg
```

#### Plot Specific Sample

```bash
quantmsioc visualize plot ibaq-distribution \
    --ibaq-path tests/examples/AE/PXD016999.1-ibaq.tsv \
    --select-column Sample_001 \
    --save-path ./plots/ibaq_sample001.svg
```

### Output

- **Format**: SVG or PDF
- **Content**: Distribution plot (histogram + kernel density estimate)
- **Axes**:
  - X-axis: log10(iBAQ intensity)
  - Y-axis: Density or frequency

### Interpretation

- **Bimodal distribution**: May indicate distinct protein abundance classes
- **Long tail**: High-abundance proteins (housekeeping, abundant structural proteins)
- **Narrow range**: Limited dynamic range, possible detection issues

### Best Practices

- Log-transform intensities for better visualization
- Compare distributions across samples to identify outliers
- Check for batch effects if distributions vary significantly

---

## kde-intensity

Plot Kernel Density Estimation (KDE) of intensity distributions across samples.

### Description

Creates overlaid KDE plots showing intensity distributions for multiple samples, enabling visual comparison of sample-to-sample variability and batch effects.

### Parameters

| Parameter        | Type    | Required | Default | Description                        |
| ---------------- | ------- | -------- | ------- | ---------------------------------- |
| `--feature-path` | Path    | Yes      | -       | Feature parquet file path          |
| `--save-path`    | Path    | Yes      | -       | Output image path (e.g., plot.svg) |
| `--num-samples`  | Integer | No       | 10      | Number of samples to plot          |

### Usage Examples

#### Plot Default Samples

```bash
quantmsioc visualize plot kde-intensity \
    --feature-path tests/examples/parquet/feature.parquet \
    --save-path ./plots/intensity_kde.svg
```

#### Plot More Samples

```bash
quantmsioc visualize plot kde-intensity \
    --feature-path ./output/feature.parquet \
    --save-path ./plots/intensity_kde_all.svg \
    --num-samples 20
```

### Output

- **Format**: SVG or PDF
- **Content**: Overlaid KDE curves for each sample
- **Axes**:
  - X-axis: log10(intensity)
  - Y-axis: Density
- **Legend**: Sample identifiers

### Interpretation

- **Overlapping curves**: Good sample-to-sample consistency
- **Shifted curves**: Potential batch effects or normalization issues
- **Different shapes**: Sample-specific technical issues

### Best Practices

- Limit to 10-20 samples for readability
- Use this plot to identify samples requiring normalization
- Compare before and after normalization

---

## peptide-distribution

Plot the distribution of peptides across proteins.

### Description

Visualizes how many peptides are identified for each protein, providing insights into protein coverage and identification confidence.

### Parameters

| Parameter        | Type    | Required | Default | Description                        |
| ---------------- | ------- | -------- | ------- | ---------------------------------- |
| `--feature-path` | Path    | Yes      | -       | Feature parquet file path          |
| `--save-path`    | Path    | Yes      | -       | Output image path (e.g., plot.svg) |
| `--num-samples`  | Integer | No       | 20      | Number of top proteins to display  |

### Usage Examples

#### Basic Example

```bash
quantmsioc visualize plot peptide-distribution \
    --feature-path tests/examples/parquet/feature.parquet \
    --save-path ./plots/peptide_per_protein.svg
```

#### Show Top 50 Proteins

```bash
quantmsioc visualize plot peptide-distribution \
    --feature-path ./output/feature.parquet \
    --save-path ./plots/peptide_per_protein_top50.svg \
    --num-samples 50
```

### Output

- **Format**: SVG or PDF
- **Content**: Bar plot showing peptide counts per protein
- **Axes**:
  - X-axis: Protein identifiers (top N by peptide count)
  - Y-axis: Number of identified peptides

### Interpretation

- **High peptide counts**: Abundant proteins with good coverage
- **Single peptide proteins**: May be less confident identifications
- **Distribution shape**: Reflects proteome complexity

### Best Practices

- Focus on top proteins for initial quality assessment
- Filter single-peptide identifications for high-confidence datasets
- Compare expected vs. observed peptide counts for key proteins

---

## box-intensity

Plot box plots of intensity distributions across samples.

### Description

Creates box plots showing the distribution of feature intensities for each sample, ideal for identifying outliers and assessing normalization quality.

### Parameters

| Parameter        | Type    | Required | Default | Description                        |
| ---------------- | ------- | -------- | ------- | ---------------------------------- |
| `--feature-path` | Path    | Yes      | -       | Feature parquet file path          |
| `--save-path`    | Path    | Yes      | -       | Output image path (e.g., plot.svg) |
| `--num-samples`  | Integer | No       | 10      | Number of samples to plot          |

### Usage Examples

#### Basic Example

```bash
quantmsioc visualize plot box-intensity \
    --feature-path tests/examples/parquet/feature.parquet \
    --save-path ./plots/intensity_boxplot.svg
```

#### Plot All Samples

```bash
quantmsioc visualize plot box-intensity \
    --feature-path ./output/feature.parquet \
    --save-path ./plots/intensity_boxplot_all.svg \
    --num-samples 50
```

### Output

- **Format**: SVG or PDF
- **Content**: Box plots for each sample
- **Axes**:
  - X-axis: Sample identifiers
  - Y-axis: log10(intensity)
- **Elements**: Box (IQR), whiskers (1.5×IQR), outliers (points)

### Interpretation

- **Aligned medians**: Good normalization
- **Similar IQR**: Consistent quantification across samples
- **Many outliers**: May indicate contamination or technical issues
- **Different ranges**: Batch effects or loading differences

### Best Practices

- Use this plot to identify samples requiring normalization
- Check for systematic differences between batches or conditions
- Compare before and after normalization to verify effectiveness
- Flag samples with unusual distributions for further investigation

---

## General Plotting Tips

### Output Formats

- **SVG**: Recommended for publications (scalable, editable)
- **PDF**: Alternative vector format (portable)
- **PNG**: Raster format (not recommended for publications)

### Color Considerations

- Plots use color-blind friendly palettes by default
- Ensure sufficient contrast for grayscale printing

### Size and Resolution

- Vector formats (SVG/PDF) scale without quality loss
- Suitable for both presentations and manuscripts

### Customization

For advanced customization beyond these commands, consider:

1. Exporting data to CSV and using custom plotting scripts
2. Using the quantmsio Python API for programmatic access
3. Importing SVG files into vector graphics editors

---

## Related Commands

- [Convert Commands](cli-convert.md) - Prepare data for visualization
- [Transform Commands](cli-transform.md) - Process data before plotting
- [Statistics Commands](cli-stats.md) - Generate numeric summaries
