#!/bin/bash

cd q_scripts
jbs=20
#---------------------------

exe=main_he

L=24
#for t1 in $(seq -10 5 60); do

t1=0
t1b=0
t1c=0

t2=0
t2b=0
t2c=0

#t3=0
for t3 in 0 10 20; do

#columnar
a1=1
a2=2

b1=0
b2=2

#for a1 in $(seq 10 14); do
#a2=$L
#b2=$L
#b1=$a1
#for b1 in $(seq 0 1 ${a1}); do

for t2 in $(seq 0 8 40); do
for t2b in $(seq 0 8 40); do
for t2c in $(seq 0 8 40); do

#t2b=$t2
#t2c=$t2

  filename=${exe}_${L}_${a1}_${t2}_${t2b}_${t2c}_${t3}.sh

  sed -e s/executable/"$exe"/g -e s/AA/"$L"/g -e s/AB/"$a1"/g -e s/AC/"$a2"/g -e s/AD/"$b1"/g -e s/AE/"$b2"/g -e s/AF/"$t1"/g -e s/AG/"$t2"/g -e s/AH/"$t3"/g -e s/AI/"$t1b"/g -e s/AJ/"$t2b"/g -e s/AK/"$t1c"/g -e s/AL/"$t2c"/g ../qsub_sample2.sh > ${filename}

  chmod 750 ${filename}

  qsub ${filename}
#  echo q_scripts/${filename}

  jbs=`expr $jbs + 1`
done
done
done
done

#----------------------------------------
cd ..
echo We inserted $jbs jobs.

