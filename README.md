# LAI_neolithic

Code for:
“Local ancestry inference identifies robust evidence of selection in Neolithic Europe”
Mies & Mathieson


Overview

This repository contains code to:

- Run local ancestry inference (LAI) using multiple methods
- Perform downstream selection analyses
- Generate all main and supplemental figures in the manuscript
- Conduct replication and tract-based analyses

Processed data required to reproduce all figures and key results are available on Zenodo:
https://doi.org/10.5281/zenodo.19684580


Quick start (figure reproduction)

Clone repository:
git clone https://github.com/gmies/LAI_neolithic.git
cd LAI_neolithic

Download Zenodo data:
https://doi.org/10.5281/zenodo.19684580

Place data as follows:
- final_figures_data      → ../final_figures_data
- final_sup_figures_data  → ../final_sup_figures_data

Uncompress data:
gunzip final_figures_data.gz
gunzip final_sup_figures_data.gz

Generate figures and tables:
Rscript figures/LAI_figures.Rmd


Repository structure

figures:
  Scripts to generate all manuscript figures
  - LAI_figures.Rmd

run_LAI:
  Pipelines for running local ancestry inference on discovery data
  - gen35_full_LAI_script.sh (7v48)
  - posterior filtering scripts
  - configurations for source sample sizes (7v1, 7v3, 7v7)

global_inference:
  Scripts for global ancestry inference
  - ADMIXTURE
  - qpAdm

LAI_analyses:
  Downstream analyses of LAI output
  - rm_LD_sliding_bins_cov_zscore.R
  - individual_ancestry_correlation.R
  - bins_correlation_matrix_names.R
  - run_rfmix_format.sh

scripts:
  Utility scripts used across workflows
  - filter_ancestry_hmm_format.R
  - rfmix_global.R
  - rfmix_output_to_mean.R

tract_analysis:
  Scripts for tract-based analyses
  - TRACTS
  - TRACTOR
  - theoretical tract length distributions

X_chromosome:
  X chromosome–specific analyses
  - LAI pipelines
  - tract length analysis (X_tractlengths.sh)
  - qq plots and replication scripts
  - run_admixture.sh
  - run_qpadm.sh
  - x_replication.sh

replication:
  Replication analyses of top signals
  - replication_1/ (ancestry_hmm; pseudohaploid data)
  - replication_2/ (RFMix; imputed data)
  - replication.R


Data availability

Processed data for all figures and analyses are available on Zenodo:
https://doi.org/10.5281/zenodo.19684580

Raw genotype data are not included in this repository.


Software and dependencies

R (>= 4.2)

Required packages:
- ggplot2
- dplyr
- cowplot
- gridExtra

Python (>= 3.9)

- numpy
- pandas
- matplotlib


External tools

The following software is required to reproduce full analyses:

- RFMix
- ADMIXTURE
- qpAdm (AdmixTools)
- ancestry_hmm
- TRACTS
- TRACTOR
- simpLAI
- Mosaic
- Gnomix
- Recomb-Mix

These must be installed separately and available in your environment.


Citation

If you use this code or data, please cite:

Mies & Mathieson, 
Local ancestry inference identifies robust evidence of selection in Neolithic Europe


Contact

For questions or issues, please open a GitHub issue or contact the authors.
