#!/bin/bash
#SBATCH --qos=debug
#SBATCH --time=60
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --cpus-per-task=2
#SBATCH --constraint=haswell
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=pcaetano@ifi.unicamp.br
#SBATCH --account=desi
#SBATCH --job-name=lightcones_postprocess

. $HOME/.bashrc
mydesienv 19.12
. envvars

conda activate lightcone 

export DIR_OUT=/global/cscratch1/sd/prc/data/mock_challenge/lightcones/v2/

srun --cpu-bind=cores multi_box2desi_cutsky_nz.py --dir_out ${DIR_OUT} config.ini

