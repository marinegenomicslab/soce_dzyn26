#!/bin/bash
# run panel design with different
# options
# written for execution on Ceres
#SBATCH --cpus-per-task=1  # ask for 1 cpu
#SBATCH --mem=10G # Maximum amount of memory this job will be given
#SBATCH --time=47:00:00 # ask that the job be allowed to run for 
#SBATCH --output=panelDesign.out # tell it where to store the output console text
#SBATCH -p short # partition to request

echo "My SLURM_JOB_ID: " $SLURM_JOB_ID

module load r/4.3.0 gcc primer3 bowtie2


# prefix
# 0 or 1 for filterNoRecomb
# 0 or 1 for filterAdj
# 75 or 150 for readlength

echo "0 0 75"
Rscript panel_design.R noRecom_noAdj_75 0 0 75
echo "0 0 150"
Rscript panel_design.R noRecom_noAdj_150 0 0 150
echo "0 1 75"
Rscript panel_design.R noRecom_Adj_75 0 1 75
echo "0 1 150"
Rscript panel_design.R noRecom_Adj_150 0 1 150
echo "1 0 75"
Rscript panel_design.R Recom_noAdj_75 1 0 75
echo "1 0 150"
Rscript panel_design.R Recom_noAdj_150 1 0 150
echo "1 1 75"
Rscript panel_design.R Recom_Adj_75 1 1 75
echo "1 1 150"
Rscript panel_design.R Recom_Adj_150 1 1 150
