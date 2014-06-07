#!/bin/bash

timestamp=`date +%y%m%d%H%M%S`
LOG=$HOME/vmc_general/log/${MACHINE}_${timestamp}

nice -n 19 $HOME/vmc_general/bin/$1 &> ${LOG} &

echo $!
