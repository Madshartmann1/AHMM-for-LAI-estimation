#!/bin/bash


#filenames and paths
VCF_FILE="l10k_n100.vcf.gz"
AHMM_INPUT="l10k_n100.panel"

# Run msprime_sim_v1.py and pipe the output to bgzip and save as sim_v1_vcf.gz
nice -19 python3 msprime_admix_sim.py /home/projects/MAAG/msprime_deme/demes/OOA_KAR_extension_admix.yaml mutation_ancestry_positions.tsv | bgzip -c > $VCF_FILE

#normalize the vcf file
/home/ctools/bcftools-1.13/bcftools norm -m -any $VCF_FILE | bgzip > norm_$VCF_FILE

# Define the number of samples and the number of populations
num_samples=500
num_populations=5

# Check if sample2pop.txt exists and the line count matches num_samples
if [ ! -f sample2pop.txt ] || [ $(wc -l < sample2pop.txt) -ne $num_samples ]; then
  # Generate the mapping file
  for i in $(seq 0 $((num_samples-1))); do
    if [ $((i / (num_samples / num_populations))) -eq 4 ]; then
      population="admixed"
    else
      population=$((i / (num_samples / num_populations)))
    fi
    printf "tsk_%d\t%s\n" "$i" "$population" >> sample2pop.txt
  done
fi

#converting into AHMM input format
python3 vcf2ahmm_fixed.py -v <(zcat norm_$VCF_FILE) -s sample2pop.txt -m 0 -g 1 > $AHMM_INPUT


#running AHMM
nice -19 ~/tools/Ancestry_HMM/src/ancestry_hmm -i $AHMM_INPUT -s ahmm.ploidy -a 4 0.17 0.33 0.5 0. -p 0 1 0.25 -p 1 1 0.25 -p 2 1 0.25 -p 3 1 0.25

#Relocating the output from AHMM
mv *.posterior posterior_files/


#Getting intervals for performance

#Theoretical intervals
nice -19 /home/ctools/R-4.1.2/bin/Rscript msprime_mutation_ancestery_intervals.R mutation_ancestry_positions.tsv theoretical_intervals.tsv
cut -f1 mutation_ancestry_positions.tsv | tail -n +2 > mutations.list


#AHMM predicted intervals
/home/ctools/R-4.1.2/bin/Rscript ahmm_mutation_estimation_intervals.R posterior_files/*

#Evaluation of performance
/home/ctools/R-4.1.2/bin/Rscript performance_3.R interval_predictions/*

