#!/bin/bash

echo submitting job to $1
ssh $1 "nice -n 19 vmc_general/local_submit.sh \"$2\""

