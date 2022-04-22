#!/bin/bash -l
#SBATCH --job-name=Discus_Opt_Cosine
#SBATCH --account=def-brami # adjust this to match the accounting group you are using to submit jobs
#SBATCH --time=6-23:00:00         # adjust this to match the walltime of your job
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32      # adjust this if you are using parallel commands
#SBATCH --mem=150000             # adjust this according to your the memory requirement per node you need
#SBATCH --mail-user=michael.melville@ryerson.ca # adjust this to match your email address
#SBATCH --mail-type=ALL

# Choose a version of MATLAB by loading a module:
module load nixpkgs/16.09
module load StdEnv/2020
module load matlab/2020b
# Remove -singleCompThread below if you are using parallel commands:
srun matlab -nodisplay -r "MAIN_DISCUS_OPT_COSINE"