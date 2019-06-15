#!/bin/bash

## job-name is just how you want the job labelled
#SBATCH --job-name=CHPD 
## output specifies the output file name
#SBATCH --output=CHPD.log
#SBATCH -p debug
## number of nodes
#SBATCH -N 1
##SBATCH --exclusive
## Can change partition from general to west
#SBATCH --partition=general
#SBATCH --mem=16GbOB
#SBATCH --time=23:59:59
#SBATCH --error=error.log


source activate rmg_env
python  /home/qin.she/Code/RMG-Py/rmg.py input.py
