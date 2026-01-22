#!/bin/bash

#to run qpadm per individual on allentoft imputed data

#cp the files here before running:
cp /project/mathilab/gmies/neolithic_selection/qpadm/0100825_allentoft_qpadm/output_data/qpadm_input* .

#rename pop for each individual:
awk '{ if ($3 == "POP") $3 = $1; printf "%-15s %-2s %-5s\n", $1, $2, $3 }' qpadm_input.ind > qpadm_input_modified.ind && mv qpadm_input_modified.ind qpadm_input.ind

module load R/4.2

Rscript ../scripts/run_sep_qpadm.R /project/mathilab/gmies/neolithic_selection/qpadm/0100825_allentoft_qpadm/output_data/qpadm_input.ind
