#tract lengths:

#rfmix:
cd ../rfmix
python  /project/mathilab/gmies/neolithic_selection/X_chromosome/LAI/scripts/analyze_rfmix_tracts.py \
--hmm-file /project/mathilab/gmies/neolithic_selection/X_chromosome/LAI/ancestry_hmm_7v7/hap/output_hmm_run_file \
--rfmix-pattern "output_rfmix_\${chr}_final.txt.gz" \
  --output-prefix rfmix_chr23 \
  --chr-start 23 \
  --chr-end 23


#simplai:
cd ../simplai
python  /project/mathilab/gmies/neolithic_selection/X_chromosome/LAI/scripts/analyze_rfmix_tracts.py \
--hmm-file /project/mathilab/gmies/neolithic_selection/X_chromosome/LAI/ancestry_hmm_7v7/hap/output_hmm_run_file \
--rfmix-pattern "output_simplai_\${chr}_final.txt.gz" \
  --output-prefix simplai_chr23 \
  --chr-start 23 \
  --chr-end 23


#ancestry hmm: (different input name)
cd ../ancestry_hmm
python  /project/mathilab/gmies/neolithic_selection/X_chromosome/LAI/scripts/analyze_rfmix_tracts.py \
--hmm-file /project/mathilab/gmies/neolithic_selection/X_chromosome/LAI/ancestry_hmm_7v7/hap/output_hmm_run_file \
--rfmix-pattern "output_ancestryhmm_\${chr}.txt.gz" \
  --output-prefix ancestryhmm_chr23 \
  --chr-start 23 \
  --chr-end 23


#mosaic:
cd ../mosaic
python  /project/mathilab/gmies/neolithic_selection/X_chromosome/LAI/scripts/analyze_rfmix_tracts.py \
--hmm-file /project/mathilab/gmies/neolithic_selection/X_chromosome/LAI/ancestry_hmm_7v7/hap/output_hmm_run_file \
--rfmix-pattern "output_mosaic_\${chr}_final.txt.gz" \
  --output-prefix mosaic_chr23 \
  --chr-start 23 \
  --chr-end 23


#plot hg tracts together from each method 
Rscript add_theoretical_tracts_5.R 4 hg ../ancestry_hmm_7v7/hap/ancestryhmm_chr23_tract_lengths_data.txt ../simplai/simplai_chr23_tract_lengths_data.txt ../rfmix/rfmix_chr23_tract_lengths_data.txt ../mosaic_7v7/mosaic_chr23_tract_lengths_data.txt ancestry_hmm simplai rfmix mosaic 

#plot farmer tracts together from each method
Rscript add_theoretical_tracts_5.R 4 farmer ../ancestry_hmm_7v7/hap/ancestryhmm_chr23_tract_lengths_data.txt ../simplai/simplai_chr23_tract_lengths_data.txt ../rfmix/rfmix_chr23_tract_lengths_data.txt ../mosaic_7v7/mosaic_chr23_tract_lengths_data.txt ancestry_hmm simplai rfmix mosaic



#by sex:

python analyze_rfmix_tracts_by_sex.py \
--hmm-file /project/mathilab/gmies/neolithic_selection/X_chromosome/LAI/ancestry_hmm_7v7/hap/output_hmm_run_file \
--rfmix-pattern "../rfmix/output_rfmix_\${chr}_final.txt.gz" \
--sex-ploidy sex_ploidy.txt \
  --output-prefix rfmix_chr23_sex \
  --chr-start 23 \
  --chr-end 23

python analyze_rfmix_tracts_by_sex.py \
--hmm-file /project/mathilab/gmies/neolithic_selection/X_chromosome/LAI/ancestry_hmm_7v7/hap/output_hmm_run_file \
--rfmix-pattern "../simplai/output_simplai_\${chr}_final.txt.gz" \
--sex-ploidy sex_ploidy.txt \
  --output-prefix simplai_chr23_sex \
  --chr-start 23 \
  --chr-end 23

python analyze_rfmix_tracts_by_sex.py \
--hmm-file /project/mathilab/gmies/neolithic_selection/X_chromosome/LAI/ancestry_hmm_7v7/hap/output_hmm_run_file \
--rfmix-pattern "../mosaic_7v7/output_mosaic_\${chr}_final.txt.gz" \
--sex-ploidy sex_ploidy.txt \
  --output-prefix mosaic_chr23_sex \
  --chr-start 23 \
  --chr-end 23
  
python analyze_rfmix_tracts_by_sex.py \
--hmm-file /project/mathilab/gmies/neolithic_selection/X_chromosome/LAI/ancestry_hmm_7v7/hap/output_hmm_run_file \
--rfmix-pattern "../ancestry_hmm_7v7/hap/output_ancestryhmm_\${chr}.txt.gz" \
--sex-ploidy sex_ploidy.txt \
  --output-prefix ancestryhmm_chr23_sex \
  --chr-start 23 \
  --chr-end 23

