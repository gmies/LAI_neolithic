#!/bin/bash

#7/7/26

module purge

module load miniconda
eval "$(/appl/miniconda/bin/conda shell.bash hook)"
conda activate gnomix

which python
python --version
python -c "import sys; print(sys.executable)"
python -m pip show numpy
python -m pip show uncertainty-calibration

######################################################
# GNOMIX
######################################################

#module load python/3.13

#mkdir gnomix
cd gnomix


chr=X

echo "Running chromosome ${chr}"

python \
/project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/11_12_25_run_samplesizes/gnomix/gnomix.py \
../recombmix/admixed.chrX.phased.vcf.gz \
chr23_output \
X \
True \
../shared_LAI_inputs/interpolated_maps/gnomix_map.chrX.tsv \
../recombmix/reference.chrX.phased.vcf.gz \
../shared_LAI_inputs/gnomix.smap \
/project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/11_12_25_run_samplesizes/gnomix/config.yaml


cd ..

