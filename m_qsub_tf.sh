#!/bin/bash

cd q_scripts
jbs=0
#---------------------------

pi=`echo "4*a(1)" | bc -l`

#ex=ampr10
#ex=swave_sq
#ex=dwave_tr2
#ex=dwave_tr2_a
#ex=dwave
ex=swave
#ex=did2
#for ex in dwave pwave fwave pip did; do
exe=main_tf_$ex
#exe=main_u1hyb2
#exe=main_he

L=24
#for L in 14 16 18 20 22 24 26 28 30; do
#for L in 14 16 18 20; do
#for L in 22 24 26 28 30; do
dd0=1

# scan two-flavor pairing

#for ap in 0 1; do
apx=1
apy=0

for p in $(seq 1 49); do

  dd0=`echo "100*s($p*$pi/100.)/c($p*$pi/100.)" | bc -l`
  dd=`echo "($dd0+.5)/1" | bc`

#for dd in $(seq 5 20 400); do
#for dd in $(seq 100 50 600); do
#for dd in $(seq 15 5 95); do
#let dd=8*10**$ddd
#echo $dd
#for dd in $(seq 330 30 600); do

  #dd=100
  t1=100
  r=3
#for r in 0 1; do

#for lth in $(seq 0 20 200); do

#for lr1 in $(seq 0 20 200); do
#lr2=$lr1
#for lr2 in $(seq `expr ${lr1} - 10` 10 `expr ${lr1} + 10`); do
#for lth in 50 150 200 250 300; do

#for phi1 in $(seq 0 20 140); do
#for phi2 in $(seq -100 40 ${phi1}); do

phi1=100
phi2=0
lth=50 #mu
lr1=20 #nk
lr2=30

#for t1b in $(seq 0 20 80); do
t1b=100
t1c=100

  filename=${exe}_${L}_${apx}_${t1b}_${dd}.sh

  sed -e s/executable/"$exe"/g -e s/AA/"$L"/g -e s/AB/"$apx"/g -e s/AC/"$apy"/g -e s/AD/"$t1"/g -e s/AE/"$dd"/g -e s/AF/"$phi1"/g -e s/AG/"$phi2"/g -e s/AH/"$lth"/g -e s/AI/"$t1b"/g -e s/AJ/"$t1c"/g -e s/AK/"$lr1"/g -e s/AL/"$lr2"/g ../qsub_sample2.sh > ${filename}

  chmod 750 ${filename}

  qsub ${filename}
#  echo q_scripts/${filename}

  jbs=`expr $jbs + 1`
done
#done
#done
#done
#done

#----------------------------------------
cd ..
echo We inserted $jbs jobs.

