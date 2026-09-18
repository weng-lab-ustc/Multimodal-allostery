# Multimodal allostery in a single-domain protein

R package and analysis workflow for studying multimodal allosteric regulation in KRAS.

This repository contains the analysis scripts, R package, and supporting data used to reproduce the analyses and figures described in:

**Multimodal allostery in a single-domain protein**

## Overview

We performed deep mutational scanning of KRAS against multiple interaction partners to characterize the energetic landscape of protein binding and allosteric regulation.

This repository contains:

* Analysis scripts for generating energy data and reproducing the manuscript figures
* The `multimodalallostery` R package containing reusable analysis and plotting functions
* Fitness scores and inferred binding free-energy changes
* Supplementary data and files required for the analyses
* Supplementary tables
* Structural movie supplementary files

## Data and figure files

Due to file-size limitations, some large files are hosted outside GitHub.

* **Supplementary data and supporting files:** [Supplementary_data](https://github.com/weng-lab-ustc/Multimodal-allostery/tree/main/Supplementary_data)
* **Scripts and processed data for reproducing manuscript figures:** [Google Drive](https://drive.google.com/drive/folders/15Yo7jXR0REYTcOFQUstYPGnHiNBJPwVX)
* **Large intermediate files:** [Zenodo](https://zenodo.org/records/21964541)

## Required software

### R package

The `multimodalallostery` package can be installed directly from GitHub:

```r
options(timeout = 600)

remotes::install_github(
  "weng-lab-ustc/Multimodal-allostery",
  dependencies = TRUE
)

library(multimodalallostery)
```

### DiMSum v1.2.9

[DiMSum](https://github.com/lehner-lab/DiMSum) was used for preprocessing deep mutational scanning sequencing data, from FASTQ files to variant-level fitness measurements.

The corresponding DiMSum scripts are provided in the repository.

### MoCHI

[MoCHI](https://github.com/lehner-lab/MoCHI) was used to fit mechanistic models to deep mutational scanning data and infer binding and folding free-energy changes from fitness measurements.

The corresponding MoCHI scripts are provided in the repository.

## Data availability

Raw sequencing data are available through the European Nucleotide Archive (ENA) under accession **PRJEB123262**.

Processed mutation-level binding free-energy datasets required for reproducing the manuscript figures are available through the [Google Drive folder](https://drive.google.com/drive/folders/15Yo7jXR0REYTcOFQUstYPGnHiNBJPwVX).

Large intermediate files are available through [Zenodo](https://zenodo.org/records/21964541).

## Figure reproduction

The scripts used to reproduce the main and supplementary figures are available in the [Google Drive folder](https://drive.google.com/drive/folders/15Yo7jXR0REYTcOFQUstYPGnHiNBJPwVX).

The figure-generation workflow is organized by manuscript figure and panel. The `multimodalallostery` R package provides reusable functions used by the figure-generation scripts.

## R package

### `multimodalallostery`

This package provides reusable functions for the analysis and visualization of the multimodal allostery workflow.

The package contains functions for data processing, statistical analysis, structural-region analysis, and figure/panel generation.

The orchestration of the complete figure workflow is handled by `run_all_figures.R`, which is provided separately from the package.

## Citation

If you use the code, package, or datasets from this repository, please cite the associated publication:

**Multimodal allostery in a single-domain protein**

[Publication / DOI to be added]

## Contact

For questions regarding the analysis or code, please open an issue in this repository or contact the authors.

