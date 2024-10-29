# dbCAN-Sum

**Script to aggregate multiple result files from the dbCAN3 web server's automated Carbohydrate-active enzyme ANnotation for my project.**

## Overview

This tool processes and aggregates dbCAN3 annotation files for carbohydrate-active enzymes (CAZymes) across multiple species or samples.

## Input Files

- **Source**: Test input files are generated from the dbCAN3 web server ([link here](https://bcb.unl.edu/dbCAN2/index.php)).
- **Preparation**: Each input file is the `overview.txt` output file, renamed to correspond to the species or sample name (e.g., `species1.txt`, `species2.txt`).
- **Configuration**: dbCAN3 configured for testing with:
  - **Input**: The test species' total protein file in `.faa` format.
  - **Databases**: Selected `HMMER_dbCAN`, `HMMER_dbCAN_sub`, and `DIAMOND` databases for the prediction.

## Aggregation Rules

The summarization process applies the following rules:

### Tool Filtering
- Only genes with results from at least **two tools** (from `HMMER_dbCAN`, `HMMER_dbCAN_sub`, `DIAMOND`) are considered.

### Annotation Extraction
- If a gene has multiple annotations, the annotation with the **highest occurrence** across the three tools is recorded.
  - Example: `GT54(82-362)+CBM94(374-526)`, `GT54_e0`, and `CBM94+GT54` would be recorded as `GT54` since it appears most frequently.
- In cases where multiple annotations have the same occurrence, the first annotation appearing in **HMMER_dbCAN_sub** is used as the record.
- **Note**: This rule may differ from certain workflows, which may recommend recording all CAZyme predictions as separate entries.

### Summarization
- Counts for major categories and families are aggregated per species or sample, using the renamed `.txt` files as identifiers.

## Usage

To run the script, use the following command:

```bash
python dbCAN_statistics_multi_en.py -in '*.txt'
