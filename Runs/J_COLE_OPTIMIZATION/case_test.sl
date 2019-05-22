#!/bin/bash -l
#SBATCH --job-name=case_test
#SBATCH --account=def-brami # adjust this to match the accounting group you are using to submit jobs
#SBATCH --time=0-02:15         # adjust this to match the walltime of your job
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1      # adjust this if you are using parallel commands
#SBATCH --mem=20000             # adjust this according to your the memory requirement per node you need
#SBATCH --mail-user=tkrebs@ryerson.ca # adjust this to match your email address
#SBATCH --mail-type=ALL

# Choose a version of MATLAB by loading a module:
module load matlab/2018a
# Remove -singleCompThread below if you are using parallel commands:
srun matlab -nodisplay -r "MAIN_test"