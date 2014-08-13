#!/bin/bash

#SBATCH --job-name=j3:140:-20
#SBATCH --output=/users/invites/sbieri/log/j3:140:-20.out
#SBATCH --error=/users/invites/sbieri/log/j3:140:-20.err
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=100
#SBATCH --partition=normal
#SBATCH --exclude=pyro,sarabi,elektra
#SBATCH --export=LD_LIBRARY_PATH=/users/invites/sbieri/lib/

srun nice -n 19 $HOME/vmc_general/bin/main_he_kag2 24 3 140 -20 200

