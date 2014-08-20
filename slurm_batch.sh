#!/bin/bash

#SBATCH --job-name=c0-10:6_-84_-10_300_0_300
#SBATCH --output=/users/invites/sbieri/log/c0-10:6_-84_-10_300_0_300.out
#SBATCH --error=/users/invites/sbieri/log/c0-10:6_-84_-10_300_0_300.err
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=100
#SBATCH --partition=normal
#SBATCH --exclude=pyro,sarabi,elektra
#SBATCH --export=LD_LIBRARY_PATH=/users/invites/sbieri/lib/

srun nice -n 19 $HOME/vmc_general/bin/main_u1compl 12 1 0 1 0 0 1 6 -84 -10 300 0 300 40 40

