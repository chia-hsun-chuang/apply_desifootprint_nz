#!/bin/bash
#SBATCH --qos=debug
#SBATCH --time=30
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --constraint=haswell
#SBATCH --job-name=test_lightcone_postprocess_single

. $HOME/.bashrc
mydesienv 19.12
. envvars

conda activate lightcone 

srun -n 1 multi_box2desi_cutsky_nz.py test.config --dir_out $DIR_OUT

