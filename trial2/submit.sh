#!/bin/bash

## job-name is just how you want the job labelled
#SBATCH --job-name=trial_stage
## output specifies the output file name
#SBATCH --output=trial2.log
## number of nodes
#SBATCH -N 2
##SBATCH --exclusive
## Can change partition from general to west
#SBATCH --partition=general
#SBATCH --mem=50Gb
##SBATCH --time=23:59:59
#SBATCH --error=error.log
#SBATCH --exclude=c5003

source activate rmg_env
python  /home/qin.she/Code/RMG-Py/rmg.py -p input.py
