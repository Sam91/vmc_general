#!/bin/bash

L=12
jbs=0
t1=100

# scan the 3-sublattice order parameter for xyz ordering
cd q_scripts

#for h in 40 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79; do
#for h in -35 -30 -25 -20 -19 -18 -17 -16 -15 -14 -13 -12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 1 2 3 4 5 6 7 8 9; do
##for h in -19 -18 -17 -16 -15 -14 -13 -12 -11 11 12 13 14 15 16 17 18 19; do
for h in $(seq -2 -1 -80); do
#for h in 1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 21 22 23 24 25 26 27 28 29 31 32 33 34 35 36 37 38 39 41 42 43 44 45 46 47 48 49; do
#for h in 110 120 130 140 150 160 170 180 190 200 220; do
#h=-${h}
#h=-40

#for eta in 5 10 15 20 25 30 35 40 45; do
#for eta in 21 22 23 24 25 26 27 28 29 30 31 32 33 34; do
#for eta in 36 37 38 39 41 42 43 44; do
#for phi1 in $(seq 0 10 90); do
#for the in $(seq 0 25 100); do
#for psi in $(seq 0 25 100); do
phi1=75
the=75
psi=50

for eta in $(seq 20 10 50); do
#for phi in $(seq 0 25 100); do
#eta=10
phi=75
#for t2 in -110 -120 -130; do
#for t2 in -10 -20 -40 -60 -80 -100 -120 -140; do
t2=0

#for r3 in -80 -60 -40 -20 20 40 60 80; do

#for Nz in 44 48 52 56 60 64 68 72; do
#for Nz in 40 44 48; do
Nz=48

#t2=0

#if [ $h == 0 ]; then
#  continue
#fi

  filename=q_${Nz}_${t1}_${t2}_${h}_${phi1}_${the}_${psi}_${eta}_${phi}.sh

  sed -e s/AA/"$L"/g -e s/AB/"$Nz"/g -e s/AC/"$t1"/g -e s/AD/"$t2"/g -e s/AE/"$h"/g -e s/AF/"$phi1"/g -e s/AG/"$the"/g -e s/AH/"$psi"/g -e s/AI/"$eta"/g -e s/AJ/"$phi"/g ../qsub_sample2.sh > ${filename}

  chmod 750 ${filename}

  qsub ${filename}
#  echo q_scripts/${filename}

  jbs=`expr $jbs + 1`
done
done
#done
#done
#done
#done

cd ..

echo We inserted $jbs jobs.

