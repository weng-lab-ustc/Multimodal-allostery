# Multimodal allostery in a single-domain protein

This repository contains the R package, analysis scripts, and supplementary data used in the study **“Multimodal allostery in a single-domain protein”**.

The project investigates the allosteric regulation of KRAS using deep mutational scanning (DMS), quantitative measurements of protein–protein interactions, and computational modeling.

## 1. Overview

This repository provides resources to reproduce the analyses and figures presented in the manuscript.

The main components include:

- An R package, `multimodalallostery`, for analyzing and visualizing allosteric effects.
- Scripts for generating manuscript figures.
- Supplementary datasets used in the analyses.
- Instructions and links for accessing additional input data and intermediate files.

The computational workflow combines DMS data processing, mechanistic modeling, and quantitative analysis of mutation-induced changes in protein folding and binding energies.

## 2. Repository structure

The main repository is organized as follows:

```text
Multimodal-allostery/
├── R/
│   └── R package source files
├── man/
│   └── Package documentation
├── Supplementary_data/
│   └── Supplementary datasets
├── DESCRIPTION
├── LICENSE
├── NAMESPACE
├── README.md
└── run_all_figure.R
```

Additional scripts and processed data required for reproducing the manuscript figures are available through the external data repositories described below.

## 3. Data and figure files

Due to the size of some input and intermediate datasets, not all files are stored directly in this GitHub repository.

The following resources are provided:

| Resource | Description | Location |
|---|---|---|
| GitHub repository | R package, source code, and supplementary data | [GitHub](https://github.com/weng-lab-ustc/Multimodal-allostery) |
| Figure scripts and processed data | Scripts and processed datasets used for manuscript figures | [Google Drive](https://drive.google.com/drive/folders/15Yo7jXR0REYTcOFQUstYPGnHiNBJPwVX) |
| F2a input data | K13- and K19-associated input files required for Figure 2a | [Google Drive](https://drive.google.com/drive/folders/1jb5x-AuQng7ncpXCuwVc-eIMF2TDpxBC) |
| Large intermediate datasets | Large-scale intermediate files | [Zenodo](https://zenodo.org/records/21964541) |
| Raw sequencing data | Raw DMS sequencing reads | ENA accession: PRJEB123262 |

The F2a input data folder contains the following files:

- `weights_Binding_K13.txt`
- `weights_Binding_K19.txt`

These files are used to generate the K13 and K19 binding-interface analyses in Figure 2a.

## 4. System requirements

### Software

The analysis was tested in the following environment:

| Component | Version / Environment |
|---|---|
| R | 4.6.1 |
| Operating systems | Windows and macOS |
| DiMSum | 1.2.9 |
| MoCHI | 0.9 |

The package and figure-generation scripts are written in R.

Linux has not been included in the tested environments reported here.

### Hardware

No specialized hardware is required. The analysis and model-fitting workflows can be executed using either CPU or GPU, depending on the available hardware and software configuration. GPU acceleration may improve computational efficiency for model-fitting tasks.

A standard desktop or laptop computer should be sufficient for the F2a plotting demonstration.

### Installation time

Installation typically takes approximately 2–3 minutes on a standard desktop computer, depending on the network connection and package dependencies.

## 5. Installation

The `multimodalallostery` package can be installed directly from GitHub.

Open R or RStudio and run:

```r
options(timeout = 600)

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github(
  "weng-lab-ustc/Multimodal-allostery",
  dependencies = TRUE
)
```

Load the package:

```r
library(multimodalallostery)
```

If the package loads successfully, the installation is complete.

Some namespace masking warnings involving functions such as `first()`, `between()`, and `last()` may appear because of overlapping imports from R packages. These warnings do not necessarily indicate an installation failure.

## 6. Demo: Reproducing Figure 2a

Figure 2a provides a representative example of the analysis and visualization functions implemented in this project.

The demo generates binding-interface residue-level median ΔΔG heatmaps for K13 and K19.

### Step 1. Download the input data

Download the following files from the [F2a input data folder](https://drive.google.com/drive/folders/1jb5x-AuQng7ncpXCuwVc-eIMF2TDpxBC):

- `weights_Binding_K13.txt`
- `weights_Binding_K19.txt`

Place these files in the directory specified by `WEIGHTS_DIR` in `run_all_figure.R`.

### Step 2. Download the analysis scripts and processed data

Download the figure-generation scripts and associated processed data from the [Google Drive repository](https://drive.google.com/drive/folders/15Yo7jXR0REYTcOFQUstYPGnHiNBJPwVX).

Maintain the directory structure expected by the script, or update the corresponding input paths to match your local directory structure.

### Step 3. Configure the local paths

Open `run_all_figure.R` and update the working directory and input/output paths to match your local environment.

In particular, check the following path variables:

```r
WEIGHTS_DIR
PREDICTIONS_FILE
ANNOTATION_5
ANNOTATION_8
CONTACT_SHELL
FITNESS_DIR
FITNESS_MERGE_DIR
MANUSCRIPT_FIGURE_DIR
OUTPUT_DIR
```

The script currently contains a working-directory setting from the original analysis environment. This should be replaced with the appropriate local path before execution.

### Step 4. Run the analysis

Run the figure-generation script in R or RStudio:

```r
source("run_all_figure.R")
```

The script contains the analysis and figure-generation procedures for multiple manuscript panels. It executes the corresponding analysis calls in the order specified in the script.

The F2a analysis uses the K13 and K19 input files described above and does not depend on the input data required for the preceding Figure 1 panels. However, running the entire script executes the other included panel analyses as well, so their required input data should also be available when running the complete script.

### Step 5. Expected outputs

The F2a analysis generates two PDF files:

- `F2a_K13.pdf`
- `F2a_K19.pdf`

They are saved under the `Figure2` subdirectory of `OUTPUT_DIR`.

The expected runtime for the F2a analysis is approximately 1 minute under the tested environment. Actual runtime may vary depending on computer performance.

## 7. Instructions for use

The `multimodalallostery` package provides functions for analyzing and visualizing mutation-induced changes in protein binding and allosteric regulation.

After installation, load the package:

```r
library(multimodalallostery)
```

For function-specific documentation, use:

```r
help(package = "multimodalallostery")
```

Alternatively, access individual function documentation using:

```r
?function_name
```

The complete figure-generation workflow is provided in `run_all_figure.R`. Additional input files and processed datasets required by individual analyses are available through the external data repositories listed above.

## 8. Reproducing the manuscript results

To reproduce the analyses and figures presented in the manuscript:

1. Clone or download this GitHub repository.
2. Install the `multimodalallostery` R package.
3. Download the corresponding figure scripts and processed data from Google Drive.
4. Download large intermediate datasets from Zenodo when required.
5. Configure the input and output paths in the relevant scripts.
6. Run the corresponding analysis scripts in R.

The supplementary datasets directly associated with the manuscript are available in the `Supplementary_data/` directory.

## 9. Data availability

Raw sequencing data generated in this study have been deposited in the European Nucleotide Archive (ENA) under accession number **PRJEB123262**.

The data are currently under controlled/private status and are scheduled for public release on **1 January 2027**. Before public release, access for editorial assessment or peer review should be arranged through the corresponding author as appropriate.

Processed datasets, figure-generation scripts, and intermediate files are available through the GitHub, Google Drive, and Zenodo repositories listed above.

## 10. License

This project is distributed under the MIT License. See the `LICENSE` file for the full license text.

## 11. Citation

If you use this package or the associated analysis workflow, please cite the manuscript:

**Multimodal allostery in a single-domain protein.**

The full citation and DOI will be added following publication.

## 12. Contact

For questions about the code, datasets, or analysis workflow, please contact the corresponding author or the project team.

**Laboratory:** Weng Lab, School of Life Sciences, University of Science and Technology of China (USTC), Hefei, China.

**GitHub:** [Multimodal-allostery](https://github.com/weng-lab-ustc/Multimodal-allostery)
