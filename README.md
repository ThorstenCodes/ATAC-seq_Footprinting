# TOBIAS Wrapper Scripts for ATAC-seq Footprinting Analysis

## Description

This repository demonstrates how the [TOBIAS](https://github.com/loosolab/TOBIAS) framework can be integrated into straightforward, transparent Bash pipelines.  
Rather than relying on workflow managers such as **Nextflow** or **Snakemake**, these scripts show how TOBIAS modules can be directly connected in a lightweight, reproducible manner.

The repository contains **two custom Bash scripts** that streamline and extend the functionality of the TOBIAS pipeline, simplifying ATAC-seq footprinting analysis:

* **`TOBIAS.sh`**  
  Automates the core TOBIAS workflow, including:
  - `ATACorrect` (Tn5 bias correction)  
  - `FootprintScores` (footprinting scoring)  
  - `BINDetect` (TF binding analysis)  
  Supports both single-sample and differential analysis modes.

* **`PLOTracks_differential_binding.sh`**  
  A post-processing script that leverages the outputs from `TOBIAS.sh` to:  
  - Generate publication-quality visualizations with `PlotTracks`  
  - Perform differential motif binding comparisons using `bedtools`  
  Designed for in-depth analysis of specific transcription factors (TFs) and genomic regions.

These scripts provide a **lightweight and reproducible alternative** to full workflow managers, while still enabling comprehensive ATAC-seq footprinting analysis.

## Features

* **Integrated Workflow**: Combines multiple TOBIAS steps into a single, executable script.  
* **Differential Analysis**: Easily compare footprinting between control and treated/genetically modified samples.  
* **Customizable Post-processing**: Focus your analysis on specific TFs and genomic regions of interest.  
* **Visualization**: Generates `TOBIAS PlotTracks` for visualizing TF footprints, ATAC-seq signals, and gene annotations in specific loci.  
* **Extended Analysis**: Includes custom functions for differential **motif analysis** using `bedtools`, providing an added layer of biological insight.  


## Prerequisites

To run these scripts, you must have the following software installed and accessible in your system's `PATH`:

* **[TOBIAS](https://github.com/loosolab/TOBIAS)** (v0.15.0 or later recommended)
* **[bedtools](https://bedtools.readthedocs.io/en/latest/)** (v2.30.0 or later recommended)  
* Standard command-line tools: `wget`, `gunzip`, `sort`, `awk`, `sed`, `uniq`  

## Installation

1. Clone this repository to your local machine:

    ```bash
    git clone https://github.com/your_username/your_repo_name.git
    cd your_repo_name
    ```

2. Make the scripts executable:

    ```bash
    chmod +x TOBIAS.sh PLOTracks_differential_binding.sh
    ```

## Usage

All scripts support `-h` or `--help` to display available arguments and usage examples. For example:

```bash TOBIAS.sh -h ```
or
```bash PLOTracks_differential_binding.sh --help```

### Script 1: `TOBIAS.sh`

This script handles the core TOBIAS workflow for ATAC-seq footprinting. You must provide the necessary input files as arguments.

**Single-sample analysis**:  
This example runs the full pipeline for all BAM files in the `/path/to/BAMs` directory.

```bash TOBIAS.sh /path/to/BAMs peaks.bed genome.fa blacklist.bed motifs.meme results 8 "" Single```

**Differential analysis**:
This example compares Control_Sample against Treated_Sample and restricts the analysis to regions specified in regions.bed.

```bash ./TOBIAS.sh /path/to/BAMs peaks.bed genome.fa blacklist.bed motifs.meme results 8 /path/to/regions.bed Differential Control_Sample Treated_Sample```

### Script 2: PLOTracks_differential_binding.sh

This script is for advanced analysis and visualization of BINDetect outputs. It requires a specific configuration of input files.

**Example 1**: Using a sample table for PlotTracks and analyzing specific TFs
This command analyzes IRF1 and IRF2 binding using the settings from sample_table.txt.

´´´bash PLOTracks_differential_binding.sh Human "IRF1 IRF2" loci.bed sample_table.txt "" "" "" /path/to/TF_outputs "Sample1 Sample2" Control_Sample bound´´´


**Example 2**: Using a TF list file and specifying bigwigs directly
This command reads a list of TFs from tfs.txt and uses a space-separated list of bigwigs and their corresponding colors for visualization.

´´´bash PLOTracks_differential_binding.sh Mouse tfs.txt loci.bed "" "/path/to/bw/*.bw" "Sample1 Sample2" "red blue" /path/to/TF_outputs "Sample1 Sample2" "" all´´´

## Contributing

If you find a bug or have an idea for an improvement, feel free to open an issue or submit a pull request. Contributions are welcome.

## Citation

## Acknowledgment / Citation

This repository provides wrapper scripts that utilize the **TOBIAS** framework for ATAC-seq footprinting analysis.  
If you use TOBIAS in your work, please cite the original publication:

Bentsen et al., *Nature Communications*, 2020.

