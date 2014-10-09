#!/bin/bash

#SBATCH --job-name=rx-20-1.3
#SBATCH --output=/users/invites/sbieri/log/rx-20-1.3.out
#SBATCH --error=/users/invites/sbieri/log/rx-20-1.3.err
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=100
#SBATCH --partition=normal
#SBATCH --exclude=pyro,sarabi,elektra
#SBATCH --export=LD_LIBRARY_PATH=/users/invites/sbieri/lib/

srun nice -n 19 $HOME/vmc_general/bin/main_u1complex 20 1 0 1 0 1 3 2 26 72 300 120 300 30 50

