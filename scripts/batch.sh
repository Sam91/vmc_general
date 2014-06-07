#!/bin/bash

runon=vega

tot=100
i=0

#for xi3 in `seq 4 2 94`; do
for xitot in `seq 8 8 100`; do

  xi3=`expr $tot - $xitot`

#  xitot=`expr $tot + $xi3 / 4`
#  xi1end=`expr 0 - $xi3`

#for xi3 in -20 -40 -60 -80 -100
#for xi3 in 0 20 40 60 80 100
#  xi1=`expr $tot + $xi2`
#  echo $xi1

  ./submit.sh ${runon} "main_u1real_scan12 8 1 0 0 0 $xitot -1 $xi3 10"
#  ./submit.sh ${runon} "main_u1real_scan12 8 1 0 0 0 -$xi3 -2 $xitot $xi3 0 10"
  sleep .1

  let "i=i+1"
done

echo Submitted $i jobs.

