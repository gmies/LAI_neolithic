#!/bin/bash

#module unload python
module load python/3.10
#module load tcltk
module load tcltk/8.6.8 

# Set paths
BASE_DIR="/project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/tract_lengths/TRACTS"
SCRIPT_DIR="$BASE_DIR/scripts/TRACTS_plotting"
#SCRIPT_DIR="$BASE_DIR/fixed"
FANCYPLOT_SCRIPT="/project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/tract_lengths/TRACTS/tracts-master/example/fancyplotting.py"

# Define all target directories
#DIRS=("ancestry_hmm" "ancestry_paths" "mosaic" "rfmix" "simplai" \
#DIRS=("ancestry_hmm_phased" "ancestry_paths_phased" "mosaic_phased" "rfmix_phased" "simplai_phased")

DIRS=("mosaic_phased" "mosaic")

# Loop through each directory
for dir in "${DIRS[@]}"; do
    FULL_PATH="${BASE_DIR}/${dir}"
    echo "Processing directory: $FULL_PATH"

    # Check if directory exists
    if [ ! -d "$FULL_PATH" ]; then
        echo "Directory $FULL_PATH does not exist, skipping..."
        continue
    fi

    cd "$FULL_PATH" || continue

    # Unzip any .bed.gz files
    echo "Unzipping .bed.gz files..."
    gunzip -f *.bed.gz 

    # Run tracts
    echo "Running run_tracts.py..."
    python "$SCRIPT_DIR/922_run_tracts.py"

    # Rename outputs to expected format
    echo "Renaming Tracts output files..."
    mv -f one_pulse_tract_length_bins one_pulse_bins
    mv -f one_pulse_sample_tract_distribution one_pulse_dat
    mv -f one_pulse_predicted_tract_distribution one_pulse_pred

    # Run fancyplotting.py
    echo "Running fancyplotting.py..."
    #python "$FANCYPLOT_SCRIPT" --name one_pulse --population-tags HG,FARMER --plot-format png

    # Recompress .bed files
    echo "Recompressing .bed files..."
    gzip -f *.bed 

    # Return to original directory
    cd - > /dev/null
done

cd ../fixed

    mv -f one_pulse_tract_length_bins one_pulse_bins
    mv -f one_pulse_sample_tract_distribution one_pulse_dat
    mv -f one_pulse_predicted_tract_distribution one_pulse_pred

python "$FANCYPLOT_SCRIPT" --name one_pulse --population-tags HG,FARMER --plot-format png


echo "All done!"
