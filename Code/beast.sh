#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --partition=comp72
#SBATCH --ntasks=32
#SBATCH --job-name=beast
#SBATCH --mail-user=pedrohenriquepezzi@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=/scrfs/storage/ppezzi/BEAST_Agalinis/beast_XXXXX_run1or2.%j.out

#Display the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on `hostname`
echo Job started at `date +"%T %a %d %b %Y"`
echo Directory is `pwd`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

#Assign path variables
WD=/scrfs/storage/ppezzi/BEAST_Agalinis
THREADS=32
SEED=$RANDOM$RANDOM

#Enter directory
cd $WD

#Load modules
module purge

#run beast
echo "Using seed: $SEED"
../beast/bin/beast -seed $SEED -threads $THREADS BEAST_XXXXXXX_input_Sigma1_100kk.xml

# Final time stamp
echo Job finished at `date +"%T %a %d %b %Y"`