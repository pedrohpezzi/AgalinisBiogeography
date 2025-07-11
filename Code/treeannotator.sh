#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --partition=comp72
#SBATCH --ntasks=32       #number of cpus to use
#SBATCH --job-name=tree_annotator
#SBATCH --mail-user=pedrohenriquepezzi@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=/scrfs/storage/ppezzi/BEAST_Agalinis/concatenated/treeannotator.%j.out

#Display the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on `hostname`
echo Job started at `date +"%T %a %d %b %Y"`
echo Directory is `pwd`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

#Assign path variables
WD=/scrfs/storage/ppezzi/BEAST_Agalinis/concatenated
THREADS=32

#enter directory
cd $WD

#load modules
module purge
module load java

#set Java memory limit to avoid crashes
export _JAVA_OPTIONS="-Xms2g -Xmx32g"

#run beast
../beast/bin/treeannotator -height mean -burnin 0 -lowMem true -file BEAST_concatenated_runs_combined.trees concatenated_annotated.tre

# Final time stamp
echo Job finished at `date +"%T %a %d %b %Y"`
