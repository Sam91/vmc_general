#!/bin/bash

#SBATCH --job-name=<jobname>
#SBATCH --output=$HOME/log/<jobname>.out
#SBATCH --error=$HOME/log/<jobname>.err
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=10
#SBATCH --partition=normal

srun nice -n 19 /users/invites/sbieri/vmc_general/bin/<cmd>
