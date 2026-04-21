# LAI_neolithic
github code for "Local ancestry inference identifies robust selection signals in ancient Neolithic populations" (Mies & Mathieson).

Data to generate main text figures, supplemental text figures, and LAI calls are available at 10.5281/zenodo.19684580. 


Directories are as follow:

figures:
  main_text_figures_scripts.Rmd is an R Markdown script you can run if you put data downloaded from directory final_figures_data on zenodo into ../final_figures_data from where the R Markdown file is and uncompress to regenerate all of the main text figures 
  supplemental_figures.Rmd is an R Markdown script you can run if you put data downloaded from directory final_sup_figures_data on zenodo into ../final_sup_figures_data from where the R Markdown file is and uncompress to regenerate all of the supplemental text figures 
  Figure4C-E_tractsplot.py is python script to generate Figure 4 panels C-E using tracts

global_inference: scripts to run admixture and qpadm on discovery dataset

run_LAI: scripts to run LAI on discovery data
  .sh files for 7v1, 7v3, 7v7, and 7v48 source sample sizes (gen35_full_LAI_script.sh)
  scripts to run posterior filtering 

  scripts: scripts called within run_LAI directory scripts to create intermediate files to run LAI methods

LAI analyses: scripts to run analyses on LAI calls
  run_rfmix_format.sh converts calls for each LAI method to RFMix format (RFMix format calls available on zenodo) 
  rm_LD_sliding_bins_cov_zscore.R performs Z-score analysis on local ancestry averages for each method
  individual_ancestry_correlation.R runs correlation of local ancestry averages for each method across individuals 
  bins_correlation_matrix_names.R runs correlation of local ancestry averages across methods 

  scripts:
    filter_ancestry_hmm_format.R 
    rfmix_global.R converts RFMix format LAI calls to global ancestry estimates for each individual
    rfmix_output_to_mean.R converts RFMix format LAI calls to local ancestry averages

tract_analysis: scripts to run TRACTS, TRACTOR, and running tract analysis with theoretical tracts

X_chromosome: .sh scripts for running LAI for each method
  X_tractlengths.sh to run tract length analysis for X chromosome
  compre_sex.R
  qq_plot scripts
  run_admixture.sh
  run_qpadm.sh
  x_replication.sh 

replication: 
  replication_1: ancestry hmm scripts to run LAI on replication dataset 1 of psuedohaploid data
  replication_2: rfmix scripts to run LAI on replication dataset 2 with imputed data
  replication.R is a script to run replication analysis of top hits
