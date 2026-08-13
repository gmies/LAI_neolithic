#tract lengths:


#simplai:
cd ../simplai
python  /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/tract_lengths/analyze_rfmix_tracts.py \
--hmm-file /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/ancestry_hmm/output_hmm_run_file \
--rfmix-pattern "output_simplai_\${chr}.txt.gz" \
--output-prefix simplai


#ancestry hmm:
cd ../ancestry_hmm
for file in *.txt.gz; do [ -f "$file" ] && mv "$file" "${file%.gz}" && gzip "${file%.gz}"; done
python  /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/tract_lengths/analyze_rfmix_tracts.py \
--hmm-file /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/ancestry_hmm/output_hmm_run_file \
--rfmix-pattern "output_ancestryhmm_\${chr}.txt.gz" \
--output-prefix ancestry_hmm

#rfmix:
cd ../rfmix
python  /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/tract_lengths/analyze_rfmix_tracts.py \
--hmm-file /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/ancestry_hmm/output_hmm_run_file \
--rfmix-pattern "output_rfmix_\${chr}.txt.gz" \
--output-prefix rfmix


#mosaic:
cd ../mosaic
python  /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/tract_lengths/analyze_rfmix_tracts.py \
--hmm-file /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/ancestry_hmm/output_hmm_run_file \
--rfmix-pattern "output_mosaic_\${chr}.txt.gz" \
--output-prefix mosaic


#gnomix:
cd ../gnomix
cd /project/mathilab/gmies/neolithic_selection/1kg_runs/results/rfmix_format_022226/gnomix/
python  /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/tract_lengths/analyze_rfmix_tracts.py \
--hmm-file /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/ancestry_hmm/output_hmm_run_file \
--rfmix-pattern "output_gnomix_\${chr}.txt.gz" \
--output-prefix gnomix


#recombmix:
cd ../recombmix
cd /project/mathilab/gmies/neolithic_selection/1kg_runs/results/rfmix_format_022226/recombmix/
python  /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/tract_lengths/analyze_rfmix_tracts.py \
--hmm-file /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/ancestry_hmm/output_hmm_run_file \
--rfmix-pattern "output_recombmix_\${chr}.txt.gz" \
--output-prefix recombmix


#previously from /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/tract_lengths

#plot hg tracts together from each method 
Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/tract_lengths/theoretical/add_theoretical_tracts_5.R 6 hg ../ancestry_hmm/ancestry_hmm_tract_lengths_data.txt ../simplai/simplai_tract_lengths_data.txt ../rfmix/rfmix_tract_lengths_data.txt ../mosaic/mosaic_tract_lengths_data.txt ../recombmix/recombmix_tract_lengths_data.txt ../gnomix/gnomix_tract_lengths_data.txt ancestry_hmm simplai rfmix mosaic recombmix gnomix 

#plot farmer tracts together from each method
Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/tract_lengths/theoretical/add_theoretical_tracts_5.R 6 farmer ../ancestry_hmm/ancestry_hmm_tract_lengths_data.txt ../simplai/simplai_tract_lengths_data.txt ../rfmix/rfmix_tract_lengths_data.txt ../mosaic/mosaic_tract_lengths_data.txt ../recombmix/recombmix_tract_lengths_data.txt ../gnomix/gnomix_tract_lengths_data.txt ancestry_hmm simplai rfmix mosaic recombmix gnomix 
