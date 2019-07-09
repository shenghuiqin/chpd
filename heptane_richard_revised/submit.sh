#!/bin/bash

## job-name is just how you want the job labelled
#SBATCH --job-name=heptane
## output specifies the output file name
#SBATCH --output=heptane.log
## number of nodes
#SBATCH -N 1
##SBATCH --exclusive
## Can change partition from general to west
#SBATCH --partition=general
#SBATCH --mem=120Gb
## number of cores, should match -n flag in python rmg command below
#SBATCH -n 12
##SBATCH --time=23:59:59
#SBATCH --error=error.log
#SBATCH --exclude=c5003

source activate rmg_env
python  /home/qin.she/Code/RMG-Py/rmg.py -p -n 12 input.py
