#!/bin/bash
# calculating He for microhaps
#SBATCH --cpus-per-task=1  # ask for 1 cpu
#SBATCH --mem=10G # Maximum amount of memory this job will be given
#SBATCH --time=2-00:00:00 # ask that the job be allowed to run for 
#SBATCH --output=calc_mh_he.out # tell it where to store the output console text
#SBATCH --export=NONE # don't start with environment called from

module load java/11.0.2

# for each file
for y in $(ls allWhatshap_*.vcf.gz)
do
	popName=${y//allWhatshap_/}
	popName=${popName//\.vcf\.gz/}
	
	# 50bp window for 75bp reads
	java -jar /project/oyster_gs_sim/CalcHe_mh_vcf.jar -v $y -w 50
	mv vcf_he_output.txt He_50_"$popName".txt

	# 125bp window for 150bp reads
	java -jar /project/oyster_gs_sim/CalcHe_mh_vcf.jar -v $y -w 125
	mv vcf_he_output.txt He_125_"$popName".txt

done
