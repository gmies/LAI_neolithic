#!/bin/bash
#pipeline to run tractor for all LAI methods

#will go through each method and first convert to rfmix v2 format
#then run tractor on this
#then do tract lengths on the results -- need to convert back to rfmix v1 format (this is simple)
#then plot and plot all together
#and gzip the files 

mkdir rfmix
cd rfmix 

#v2 format:
python ../make_input_file.py ../../7v7/rfmix_format/rfmix/output_rfmix_

gunzip *.gz

#tractor:
for chr in {1..22}; do
python /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/unkink_2way_mspfile.py --msp rfmixv2_${chr}
awk '{for(i=7; i<=NF; i++) printf "%s%s", $i, (i<NF?" ":"\n")}' rfmixv2_${chr}.Unkinked.msp.tsv > rfmixv2_${chr}.ancestry_only.txt
done

#tract length analysis:
python /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/tract_lengths/analyze_rfmix_tracts.py \
--hmm-file /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/ancestry_hmm/output_hmm_run_file \
--rfmix-pattern "rfmixv2_\${chr}.ancestry_only.txt" \
--output-prefix phased_tracts

gzip rfmixv*

cd ..

mkdir simplai
cd simplai 

#v2 format:
python ../make_input_file.py ../../7v3/rfmix_format/simplai/output_simplai_

gunzip *.gz

#tractor:
for chr in {1..22}; do
python /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/unkink_2way_mspfile.py --msp rfmixv2_${chr}
awk '{for(i=7; i<=NF; i++) printf "%s%s", $i, (i<NF?" ":"\n")}' rfmixv2_${chr}.Unkinked.msp.tsv > rfmixv2_${chr}.ancestry_only.txt
done

#tract length analysis:
python /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/tract_lengths/analyze_rfmix_tracts.py \
--hmm-file /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/ancestry_hmm/output_hmm_run_file \
--rfmix-pattern "rfmixv2_\${chr}.ancestry_only.txt" \
--output-prefix phased_tracts

gzip rfmixv*

cd ..


mkdir mosaic
cd mosaic 

#v2 format:
python ../make_input_file.py ../../7v3/rfmix_format/mosaic/output_mosaic_

gunzip *.gz

#tractor:
for chr in {1..22}; do
python /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/unkink_2way_mspfile.py --msp rfmixv2_${chr}
awk '{for(i=7; i<=NF; i++) printf "%s%s", $i, (i<NF?" ":"\n")}' rfmixv2_${chr}.Unkinked.msp.tsv > rfmixv2_${chr}.ancestry_only.txt
done

#tract length analysis:
python /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/tract_lengths/analyze_rfmix_tracts.py \
--hmm-file /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/ancestry_hmm/output_hmm_run_file \
--rfmix-pattern "rfmixv2_\${chr}.ancestry_only.txt" \
--output-prefix phased_tracts

gzip rfmixv*

cd ..


mkdir ancestry_paths
cd ancestry_paths 

#v2 format:
python ../make_input_file.py /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/ancestry_paths/output_ancestrypaths_

gunzip *.gz

#tractor:
for chr in {1..22}; do
python /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/unkink_2way_mspfile.py --msp rfmixv2_${chr}
awk '{for(i=7; i<=NF; i++) printf "%s%s", $i, (i<NF?" ":"\n")}' rfmixv2_${chr}.Unkinked.msp.tsv > rfmixv2_${chr}.ancestry_only.txt
done

#tract length analysis:
python /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/tract_lengths/analyze_rfmix_tracts.py \
--hmm-file /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/ancestry_hmm/output_hmm_run_file \
--rfmix-pattern "rfmixv2_\${chr}.ancestry_only.txt" \
--output-prefix phased_tracts

gzip rfmixv*

cd ..


mkdir ancestry_hmm_haplotypes
cd ancestry_hmm_haplotypes 

#v2 format:
python ../make_input_file.py ../../7v7/rfmix_format/ancestry_hmm_haplotypes/output_ancestryhmm_

gunzip *.gz

#tractor:
for chr in {1..22}; do
python /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/unkink_2way_mspfile.py --msp rfmixv2_${chr}
awk '{for(i=7; i<=NF; i++) printf "%s%s", $i, (i<NF?" ":"\n")}' rfmixv2_${chr}.Unkinked.msp.tsv > rfmixv2_${chr}.ancestry_only.txt
done

#tract length analysis:
python /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/tract_lengths/analyze_rfmix_tracts.py \
--hmm-file /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/ancestry_hmm/output_hmm_run_file \
--rfmix-pattern "rfmixv2_\${chr}.ancestry_only.txt" \
--output-prefix phased_tracts

gzip rfmixv*

cd ..



mkdir gnomix
cd gnomix

#v2 format:
python ../make_input_file.py ../../7v7/rfmix_format/gnomix/output_gnomix_

gunzip *.gz

#tractor:
for chr in {1..22}; do
python /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/unkink_2way_mspfile.py --msp rfmixv2_${chr}
awk '{for(i=7; i<=NF; i++) printf "%s%s", $i, (i<NF?" ":"\n")}' rfmixv2_${chr}.Unkinked.msp.tsv > rfmixv2_${chr}.ancestry_only.txt
done

#tract length analysis:
python /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/tract_lengths/analyze_rfmix_tracts.py \
--hmm-file /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/ancestry_hmm/output_hmm_run_file \
--rfmix-pattern "rfmixv2_\${chr}.ancestry_only.txt" \
--output-prefix phased_tracts

gzip rfmixv*

cd ..

mkdir recombmix
cd recombmix

#v2 format:
python ../make_input_file.py ../../7v7/rfmix_format/recombmix/output_recombmix_

gunzip *.gz

#tractor:
for chr in {1..22}; do
python /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/unkink_2way_mspfile.py --msp rfmixv2_${chr}
awk '{for(i=7; i<=NF; i++) printf "%s%s", $i, (i<NF?" ":"\n")}' rfmixv2_${chr}.Unkinked.msp.tsv > rfmixv2_${chr}.ancestry_only.txt
done

#tract length analysis:
python /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/tract_lengths/analyze_rfmix_tracts.py \
--hmm-file /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/ancestry_hmm/output_hmm_run_file \
--rfmix-pattern "rfmixv2_\${chr}.ancestry_only.txt" \
--output-prefix phased_tracts

gzip rfmixv*

cd ..



#plot all together 

mkdir revision_plots
cd revision_plots 

module load R/4.4

#plot hg tracts together from each method 
#Rscript ../plots/tract_lengths_ancestry.R 5 hg ../ancestry_hmm_haplotypes/phased_tracts_tract_lengths_data.txt ../ancestry_paths/phased_tracts_tract_lengths_data.txt ../mosaic/phased_tracts_tract_lengths_data.txt ../simplai/phased_tracts_tract_lengths_data.txt ../rfmix/phased_tracts_tract_lengths_data.txt ancestry_hmm_haplotypes ancestry_paths mosaic simplai rfmix

#plot farmer tracts together from each method
#Rscript ../plots/tract_lengths_ancestry.R 5 farmer ../ancestry_hmm_haplotypes/phased_tracts_tract_lengths_data.txt  ../ancestry_paths/phased_tracts_tract_lengths_data.txt ../mosaic/phased_tracts_tract_lengths_data.txt ../simplai/phased_tracts_tract_lengths_data.txt ../rfmix/phased_tracts_tract_lengths_data.txt ancestry_hmm_haplotypes ancestry_paths mosaic simplai rfmix

#plot hg tracts together from each method 
Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/tract_lengths/plots/tract_lengths_ancestry.R 7 hg ../ancestry_hmm_haplotypes/phased_tracts_tract_lengths_data.txt ../ancestry_paths/phased_tracts_tract_lengths_data.txt ../mosaic/phased_tracts_tract_lengths_data.txt ../simplai/phased_tracts_tract_lengths_data.txt ../rfmix/phased_tracts_tract_lengths_data.txt ../gnomix/phased_tracts_tract_lengths_data.txt ../recombmix/phased_tracts_tract_lengths_data.txt ancestry_hmm_haplotypes ancestry_paths mosaic simplai rfmix gnomix recombmix

#plot farmer tracts together from each method
Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/tract_lengths/plots/tract_lengths_ancestry.R 7 farmer ../ancestry_hmm_haplotypes/phased_tracts_tract_lengths_data.txt  ../ancestry_paths/phased_tracts_tract_lengths_data.txt ../mosaic/phased_tracts_tract_lengths_data.txt ../simplai/phased_tracts_tract_lengths_data.txt ../rfmix/phased_tracts_tract_lengths_data.txt ../gnomix/phased_tracts_tract_lengths_data.txt ../recombmix/phased_tracts_tract_lengths_data.txt ancestry_hmm_haplotypes ancestry_paths mosaic simplai rfmix gnomix recombmix
