#!/bin/bash

#echo $1
#echo $2

sed "s/<jobname>/$1/g;s/<cmd>/$2/g" $HOME/vmc_general/scripts/slurm_skeleton.sh > slurm_batch.sh

sbatch slurm_batch.sh


