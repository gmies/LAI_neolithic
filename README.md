# LAI_neolithic

Code for:
“Local ancestry inference identifies robust selection signals in ancient Neolithic populations”
Mies & Mathieson

Overview

This repository contains code to:
1. Run local ancestry inference (LAI) using multiple methods
2. Perform downstream selection analyses
3. Generate all main and supplemental figures in the manuscript
4. Conduct replication and tract-based analyses

Processed data required to reproduce all figures and key results are available on Zenodo:
https://doi.org/10.5281/zenodo.19684580


Quick start (figure reproduction)
 1. Clone repository:
    git clone https://github.com/gmies/LAI_neolithic.git
    cd LAI_neolithic

 3. Download Zenodo data:
    https://doi.org/10.5281/zenodo.19684580

 5. Place data as follows:
    final_figures_data      → ../final_figures_data
    final_sup_figures_data  → ../final_sup_figures_data

 7. Uncompress data:
    gunzip final_figures_data.gz
    gunzip final_sup_figures_data.gz

 9. Generate figures:
     Rscript figures/main_text_figures_scripts.Rmd
    Rscript figures/supplemental_figures.Rmd


Repository structure:

figures/
  Scripts to generate all manuscript figures
  - main_text_figures_scripts.Rmd
  - supplemental_figures.Rmd
  - Figure4C-E_tractsplot.py

run_LAI/
  Pipelines for running local ancestry inference on discovery data
  - gen35_full_LAI_script.sh
  - scripts for posterior filtering
  - configurations for different sample sizes (7v1, 7v3, 7v7, 7v48)

global_inference/
  Scripts for global ancestry inference
  - ADMIXTURE
  - qpAdm

LAI_analyses/
  Downstream analyses of LAI output
  - rm_LD_sliding_bins_cov_zscore.R
  - individual_ancestry_correlation.R
  - bins_correlation_matrix_names.R
  - run_rfmix_format.sh

scripts/
  Utility scripts used across workflows
  - filter_ancestry_hmm_format.R
  - rfmix_global.R
  - rfmix_output_to_mean.R

tract_analysis/
  Scripts for tract-based analyses
  - TRACTS
  - TRACTOR
  - theoretical tract simulations

X_chromosome/
  X chromosome–specific analyses
  - LAI pipelines
  - tract length analysis (X_tractlengths.sh)
  - qq plots and replication scripts
  - run_admixture.sh
  - run_qpadm.sh
  - x_replication.sh

replication/
  Replication analyses of top signals
  - replication_1/ (ancestry_hmm; pseudohaploid data)
  - replication_2/ (RFMix; imputed data)
  - replication.R


Data availability
Processed data for all figures and analyses are available on Zenodo:
https://doi.org/10.5281/zenodo.19684580
Due to size and access restrictions, raw genotype data are not included in this repository.


Software and dependencies
R (≥ 4.2)
Required packages include (not exhaustive):
ggplot2
dplyr
cowplot
gridExtra

Python (≥ 3.9)
numpy
pandas
matplotlib


External tools
The following software is required to reproduce full analyses:

RFMix
ADMIXTURE
qpAdm (AdmixTools)
ancestry_hmm
TRACTS / TRACTOR
simpLAI
Mosaic

These must be installed separately and available in your environment.


Citation
If you use this code or data, please cite:
Mies & Mathieson (2026)
Local ancestry inference identifies robust selection signals in ancient Neolithic populations


Contact
For questions or issues, please open a GitHub issue or contact the authors.
