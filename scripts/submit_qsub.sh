#!/bin/bash

#echo $1
#echo $2

sed "s/<cmd>/$1/g" $HOME/vmc_general/scripts/qsub_skeleton.sh > qsub_batch.sh

qsub qsub_batch.sh


