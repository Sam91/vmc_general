#!/bin/bash

cd q_scripts
jbs=0
#---------------------------

ex=u1hyb2

#for ex in dwave pwave fwave pip did; do
exe=main_$ex
#exe=main_u1hyb2
#exe=main_he

L=12
#for L in 12 14 16 18 20; do
#for L in 24; do

# scan two-flavor pairing

#for ap in 0 1; do
ap=1

for h in $(seq 5 10 105); do
#let dd=8*10**$ddd
#echo $dd
#for dd in $(seq 330 30 600); do

#  dd=100
  t1=100

a1=1 #spiral angle in y-direction
a2=2
b1=0 #spiral angle in x-direction
b2=2

#for t1b in $(seq -20 20 20); do
for j1 in $(seq 1 2 13); do
for j1b in $(seq 1 2 13); do
j1c=$j1b
t1b=0
#j1b=0

  filename=${exe}_${L}_${ap}_${h}_${t1b}_${a2}_${b2}_${j1}_${j1b}.sh

  sed -e s/executable/"$exe"/g -e s/AA/"$L"/g -e s/AB/"$ap"/g -e s/AC/"$t1"/g -e s/AD/"$t1b"/g -e s/AE/"$h"/g -e s/AF/"$a1"/g -e s/AG/"$a2"/g -e s/AH/"$b1"/g -e s/AI/"$b2"/g -e s/AJ/"$j1"/g -e s/AK/"$j1b"/g -e s/AL/"$j1c"/g ../qsub_sample2.sh > ${filename}

  chmod 750 ${filename}

  qsub ${filename}
#  echo q_scripts/${filename}

  jbs=`expr $jbs + 1`
done
done
done
#done

#----------------------------------------
cd ..
echo We inserted $jbs jobs.

