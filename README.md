# Multimodal allostery in a single-domain protein

This repository contains the R package, analysis scripts, and supplementary data used in the study **“Multimodal allostery in a single-domain protein”**.

The project investigates the allosteric regulation of KRAS using deep mutational scanning (DMS), quantitative measurements of protein–protein interactions, and computational modeling. By comparing the effects of mutations on distinct binding interfaces, this study provides a framework for characterizing multimodal allosteric regulation in a single-domain protein.

## 1. Overview

This repository provides computational resources for reproducing the analyses and figures presented in the manuscript.

The resources include:

- An R package (`multimodalallostery`) containing functions for data analysis and visualization.
- Analysis scripts used to generate the manuscript figures.
- Instructions for installing the software and reproducing selected analyses.

The analysis workflow integrates deep mutational scanning data, quantitative protein–protein interaction measurements, and computational modeling to characterize the allosteric landscape of KRAS.

## 2. Repository structure

The GitHub repository is organized as follows:

```text
Multimodal-allostery/
├── R/                       # R functions for data analysis and visualization
├── man/                     # R package documentation
├── Supplementary_data/      # Supplementary scripts for DiMSum and MoCHI
├── DESCRIPTION              # R package metadata and dependencies
├── NAMESPACE                # Exported functions and imports
├── LICENSE                  # MIT License
└── README.md                # Project documentation
```

The main figure-generation script, `run_all_figure.R`.

The `multimodalallostery` R package provides reusable functions for the analyses, while the external analysis scripts contain the procedures used to process data and generate the manuscript figures.

## 3. Data and analysis resources

The data and analysis resources are distributed across the following locations.

### 3.1. GitHub repository

The R package and supplementary script are available in this repository:

[Multimodal-allostery on GitHub](https://github.com/weng-lab-ustc/Multimodal-allostery)

### 3.2. Analysis data

The analysis data used to generate the manuscript figures are available through Google Drive:

[Analysis data](https://drive.google.com/drive/folders/15Yo7jXR0REYTcOFQUstYPGnHiNBJPwVX)

### 3.3. Input files for Figure 2a

The following two files are required for the Figure 2a demonstration:

- `weights_Binding_K13.txt`
- `weights_Binding_K19.txt`

These files are available in the following Google Drive folder:

[Figure 2a input files](https://drive.google.com/drive/folders/1jb5x-AuQng7ncpXCuwVc-eIMF2TDpxBC)

### 3.4. Supplementary script

Supplementary script for DiMSum and MoCHI are included in the `Supplementary_script` directory of this repository:

[Supplementary_script](https://github.com/weng-lab-ustc/Multimodal-allostery/tree/main/Supplementary_script)

## 4. System requirements

The software and computational environment used in this study are listed below.

| Component | Version / Environment |
| :--- | :--- |
| R | [4.6.1](https://www.r-project.org/) |
| Operating systems | Windows and macOS |
| DiMSum | [1.2.9](https://github.com/lehner-lab/DiMSum/releases/tag/v1.2.9) |
| MoCHI | [0.9](https://github.com/lehner-lab/MoCHI/releases/tag/v0.9) |

### Hardware requirements

No specialized hardware is required.

The analysis and model-fitting workflows can be executed using either CPU or GPU, depending on the available hardware and software configuration. GPU acceleration may improve computational efficiency for model-fitting tasks.

A standard desktop or laptop computer is sufficient for running the Figure 2a plotting demonstration.

### Software dependencies

The R package dependencies are specified in the `DESCRIPTION` file and can be installed automatically during package installation.

Additional R packages used by the figure-generation scripts include `data.table`, `krasddpcams`, and `wlab.block`.

These packages can be installed from R using `install.packages()`.

## 5. Installation

The `multimodalallostery` R package can be installed directly from GitHub.

First, install the `remotes` package if it is not already available:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
```

Install `multimodalallostery` and its R dependencies:

```r
options(timeout = 600)

remotes::install_github(
  "weng-lab-ustc/Multimodal-allostery",
  dependencies = TRUE
)
```

Load the package:

```r
library(multimodalallostery)
```

Installation typically takes approximately 2–3 minutes on a standard desktop computer, depending on the network connection and package dependencies.

## 6. Demo: Reproducing Figure 2a

This demonstration reproduces the Figure 2a analysis using the K13 and K19 binding-energy datasets.

### 6.1. Download input files

Download the following two input files from the [Figure 2a input data folder](https://drive.google.com/drive/folders/1jb5x-AuQng7ncpXCuwVc-eIMF2TDpxBC):

- `weights_Binding_K13.txt`
- `weights_Binding_K19.txt`

Place the downloaded files in the corresponding input directory used by the analysis script.

### 6.2. Obtain the analysis script

Download `run_all_figure.R` from the [R](https://github.com/weng-lab-ustc/Multimodal-allostery/).

Configure the working directory and file paths in the script to match the local locations of the downloaded files and other required resources.

The script contains the analysis and figure-generation procedures for multiple manuscript panels. It executes the corresponding analysis calls in the order specified in the script.

The F2a analysis uses the K13 and K19 input files described above and does not depend on the input data required for the preceding Figure 1 panels. However, running the entire script executes the other included panel analyses as well. Therefore, all input data required by those analyses should also be available when running the complete script.

### 6.3. Run the analysis

After configuring the paths and making the required input files and software dependencies available, run the script in R or RStudio:

```r
source("run_all_figure.R")
```

### 6.4. Expected output

The Figure 2a analysis generates the following PDF files:

- `F2a_K13.pdf`
- `F2a_K19.pdf`

The output files are saved in the `Figure2` subdirectory of the specified output directory.

The Figure 2a analysis takes approximately 1 minute on a standard desktop computer.

## 7. Instructions for use

The `multimodalallostery` R package provides reusable functions for the analyses and visualization procedures presented in this study.

After installation, load the package:

```r
library(multimodalallostery)
```

To view the package documentation and available functions:

```r
help(package = "multimodalallostery")
```

Documentation for individual functions can be accessed using:

```r
?function_name
```

The associated datasets provided in this repository and through Google Drive can be used for reproducing the analyses presented in the manuscript.

Users interested in applying the analysis framework to their own datasets should refer to the function documentation and corresponding analysis scripts for information about input data formats and analysis procedures.

## 8. Reproducing manuscript results

The analysis scripts and datasets required to reproduce the manuscript figures and results are provided through this GitHub repository and the external Google Drive resources described above.

The following resources should be used together:

- **R package:** functions for data analysis and visualization.
- **Analysis scripts:** procedures for processing data and generating manuscript figures.
- **Analysis data:** available through Google Drive.

For reproduction of the complete manuscript analyses, download the relevant scripts and their corresponding input and intermediate datasets from the provided locations. Configure the local file paths in the analysis scripts before execution.


The Figure 2a demonstration described above provides a relatively small example for testing the installation and basic functionality of the software.

## 9. Data availability

Raw sequencing data have been deposited in the European Nucleotide Archive (ENA) under accession number **PRJEB123262**.

The data are currently under restricted access and are scheduled for public release on January 1, 2027. Access for editorial assessment can be arranged through the corresponding author.

Analysis data used in this study are available through [Google Drive](https://drive.google.com/drive/folders/15Yo7jXR0REYTcOFQUstYPGnHiNBJPwVX). 

## 10. License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## 11. Citation

If you use this software or the associated datasets, please cite the corresponding publication:

**Multimodal allostery in a single-domain protein.**

The full citation and DOI will be added following publication.

## 12. Contact

For questions about the software, analysis scripts, or datasets, please contact:

**Chenchun Weng**  
School of Life Sciences  
University of Science and Technology of China (USTC)  
Hefei, China
