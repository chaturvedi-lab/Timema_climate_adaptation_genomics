#!/bin/bash
#SBATCH -n 12
#SBATCH -N 1
#SBATCH -t 300:00:00
#SBATCH -p usubio-kp
#SBATCH -A usubio-kp
#SBATCH -J baypass

cd /uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/analyses/baypass/pcvars_run/outfiles/

#load appropriate modules
module load gcc
module load gsl
module load hdf5
module load perl

perl ./forkRunBaypass.pl 10
