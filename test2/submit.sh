#!/bin/bash

## job-name is just how you want the job labelled
#SBATCH --job-name=CHPD 
## output specifies the output file name
#SBATCH --output=CHPD.log
## number of nodes
#SBATCH -N 1
##SBATCH --exclusive
## Can change partition from general to west
#SBATCH --partition=west 
<<<<<<< HEAD
#SBATCH --mem=50Gb
=======
#SBATCH --mem=30Gb
>>>>>>> 09592290711a065c722241bb605aafe966cfe3b5
##SBATCH --time=100:59:59
#SBATCH --error=error.log
#SBATCH --exclude=c5003

source activate rmg_env
python  /home/qin.she/Code/RMG-Py/rmg.py -p input.py
