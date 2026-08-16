# Multimodal allostery in a single-domain protein

This repository contains the analysis scripts used to reproduce the figures and
analyses described in:

"Multimodal allostery in a single-domain protein"

## Overview

We performed deep mutational scanning of KRAS against multiple interaction
partners to characterize the energetic landscape of protein binding and
allosteric regulation.

This repository contains:
- Fitness scores, inferred free energy changes and required miscellaneous files in the [Supplementary_data] directory 
- Structural movie supplementary files 
- Script for reproducing figures in the manuscript


## Required software

### R package(devtools::install_github("weng-lab-ustc/Multimodal-allostery"))

### DiMSum v1.2.9 (pipeline for pre-processing deep mutational scanning data i.e. FASTQ to fitness)
### MoCHI (tool to fit mechanistic models to deep mutational scanning data i.e. fitness to free energy changes)

## Data availability

Raw sequencing data are available through ENA(PRJEB123262).

Processed mutation-level binding free energy datasets required for figure
reproduction are provided in the [Supplementary_data] directory.

Large intermediate files in https://zenodo.org/records/21964541?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6ImRkNGUwNTQyLTAyNTctNDA5OC1hOTg3LWI0MWEwYmJlYWViZiIsImRhdGEiOnt9LCJyYW5kb20iOiJiMTFmZTg0NGFiYWFlNjI4MDBhZTllNWUwMjNiMzQ3MiJ9.KcsbLvolTHIo_qLRdpRhKzeC0keRtmWQRy3p1wdxsR2fklXCUS2MNxVZvD7GvvWJhDobe2oP7AYR1x7ZCgYeQQ.

# multimodalallostery

Reusable functions consolidated from the Multimodal-allostery Figure/Panel workflows.

Panel orchestration remains outside the package in `run_all_figures.R`. Data are not bundled with the package.
