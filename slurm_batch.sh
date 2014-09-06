#!/bin/bash

#SBATCH --job-name=chain-90--1
#SBATCH --output=/users/invites/sbieri/log/chain-90--1.out
#SBATCH --error=/users/invites/sbieri/log/chain-90--1.err
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=100
#SBATCH --partition=normal
#SBATCH --exclude=pyro,sarabi,elektra
#SBATCH --export=LD_LIBRARY_PATH=/users/invites/sbieri/lib/

srun nice -n 19 $HOME/vmc_general/bin/main_u1chain 90 1 -1 100 0 0 100 100

