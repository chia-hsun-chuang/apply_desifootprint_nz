#!/bin/bash
#SBATCH --qos=debug
#SBATCH --time=30
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --cpus-per-task=2
#SBATCH --constraint=haswell
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=pcaetano@ifi.unicamp.br
#SBATCH --account=desi
#SBATCH --job-name=lightcone_v3_postprocess

. $HOME/.bashrc
mydesienv 19.12
. envvars

conda activate lightcone 

srun --cpu-bind=cores multi_box2desi_cutsky_nz.py config.ini --dir_out /global/cscratch1/sd/prc/data/mock_challenge/lightcones/v2/

