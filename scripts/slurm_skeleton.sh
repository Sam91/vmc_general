#!/bin/bash

#SBATCH --job-name=<jobname>
#SBATCH --output=/users/invites/sbieri/log/<jobname>.out
#SBATCH --error=/users/invites/sbieri/log/<jobname>.err
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=100
#SBATCH --partition=normal
#SBATCH --exclude=pyro,sarabi,elektra
#SBATCH --export=LD_LIBRARY_PATH=/users/invites/sbieri/lib/

srun nice -n 19 $HOME/vmc_general/bin/<cmd>

