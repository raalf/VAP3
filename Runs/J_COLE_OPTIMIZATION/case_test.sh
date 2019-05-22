#!/bin/bash

#PBS -N case_test
#PBS -m bea
#PBS -M jac064@bucknell.edu
#PBS -l pmem=20gb
#PBS -l pvmem=20gb
#PBS -l nodes=1:ppn=1
#PBS -j oe

cd $PBS_O_WORKDIR

# load the matlab module

module load matlab

matlab -nosplash -nodisplay -r MAIN_test