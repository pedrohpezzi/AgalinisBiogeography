#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --partition=tres72
#SBATCH --ntasks=1
#SBATCH --job-name=r8s_mrbayes_final
#SBATCH --mail-user=pedrohenriquepezzi
#SBATCH --mail-type=ALL
#SBATCH --output=/scrfs/storage/ppezzi/r8s/r8s_mrbayes_final.%j.out

#Display the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on `hostname`
echo Job started at `date +"%T %a %d %b %Y"`
echo Directory is `pwd`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

#Assign path variables
DIRECTORY=/scrfs/storage/ppezzi/r8s

#Load modules
module purge

#Enter directory
cd $DIRECTORY

#run r8s
r8s1.81/src/r8s -b -f AgalinisMrBayesRooted_FinalRun.nex

# Final time stamp
echo Job finished at `date +"%T %a %d %b %Y"`