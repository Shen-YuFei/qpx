# Quick Start

## Quick Start Flow

![Quick Start Flow](images/workflow2.png)

Get started with qpx in a few simple steps: from installation to data conversion, transformation, and visualization.

## Installation

```bash
pip install qpx
```

## Basic Usage

```bash
# View available commands
qpxc --help

# Convert MaxQuant data
qpxc convert maxquant-psm \
    --msms-file msms.txt \
    --output-folder ./output

# Transform to absolute expression
qpxc transform ae \
    --ibaq-file ibaq.tsv \
    --sdrf-file metadata.sdrf.tsv \
    --output-folder ./output

# Generate visualizations
qpxc visualize plot ibaq-distribution \
    --ibaq-path ./output/ae.parquet \
    --save-path ./plots/distribution.svg
```

[View complete CLI documentation →](cli-reference.md)

