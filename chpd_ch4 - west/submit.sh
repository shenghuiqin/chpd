#!/bin/bash

## job-name is just how you want the job labelled
#SBATCH --job-name=CHPD_CH4
## output specifies the output file name
#SBATCH --output=CHPD_CH4.log
## number of nodes
#SBATCH -N 2
##SBATCH --exclusive
## Can change partition from general to west
#SBATCH --partition=west
#SBATCH --mem=100Gb
##SBATCH --time=123:59:59
#SBATCH --error=error.log
#SBATCH --exclude=c5003

source activate rmg_env
python  /home/qin.she/Code/RMG-Py/rmg.py -p input.py
