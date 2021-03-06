#!/usr/bin/python
from subprocess import call
import pickle

filename = 'task_x_3'

task={}
i=0
kk=0
r=0

#for k2 in range(-600,601,20):
#  for k1 in range(20,601,20):

for L in range(16,21,2):
    #n=-1
#  for n in [-1, 0, 1, 2, 3]:
#for r3 in range(1,97,2):
 for ap in [0, 1]:
  for ex in [-1, 0, 1, 2, 3]:
#  for a2 in [0, 120, 200, 300]:
#for j2 in range(-60,101,10):
#  for j1 in range(-60,201,10):
#for j1 in range(0,601,10):
#  for j2 in range(-60,201,10):
    #j3=400
    #j2 = 30
    #k1=400
    #sgn=1
  #for sgn in [0,1]:

    #parameters for the A-phase
    r1 = 2
    r2 = 26
    r3 = 72
    a1 = 300
    a2 = 120
    a3 = 300

    #r1 = 98
    #r2 = 2
    #r3 = 0

#    t = "q3:"+str(j1)+":"+str(j2)
    t = "rx-"+ str(L) + "-" + str(ap) +"."+ str(ex)
#    r1 = 2
#    r2 = 100 - r1 - r3

    #cmd = "main_u1kagome_scan23_c 12 1 0 1 0 1 1 100 "+ str(k) +" 50 "+ str(sgn) +" 0 0 200 300 8"
    #cmd = "main_u1kagome_scan23_c 12 1 0 1 0 1 1 100 50 50 "+ str(sgn) +" 0 0 "+ str(5*k) +" 300 10"
    #cmd = "main_u1compl 8 1 0 1 0 0 1 78 13 0 "+ str(k1) +" "+ str(k2) +" 0 10"
    #cmd = "main_u1compl 8 1 0 1 0 0 1 77 23 0 "+ str(k1) +" 0 0 10"
    #cmd = "main_u1kagome_scan12_c 12 1 0 1 1 0 1 100 "+ str(k) +" 50 "+ str(sgn) +" 0 300 0 0 8"
    #cmd = "main_he_kag 20 0 -20 "+ str(j2) +" "+ str(j3) +" 200"
    #cmd = "main_he_kag2 24 3 "+ str(j1) +" "+ str(j2) +" 200"
    #cmd = "main_he_kag 20 2 -20 "+ str(j2) +" "+ str(j3) +" 200"
    #cmd = "main_u1chain "+ str(L) +" 1 "+ str(n) +" 100 0 0 100 100"
    cmd = "main_u1complex "+ str(L) +" "+ str(ap) +" 0 1 0 1 "+ str(ex) +" "+ str(r1) +" "+ str(r2) +" "+ str(r3) +" "+ str(a1) +" "+ str(a2) +" "+ str(a3) +" 30 50"
    #cmd = "main_u1realex "+ str(L) +" "+ str(ap) +" 0 1 0 0 "+ str(ex) +" "+ str(r1) +" "+ str(r2) +" "+ str(r3) +" 50 50"

    task[ t ] = [cmd, '', 'new', '' ] #command, host, status (1: unsubmitted, 0: submitted), pid

for t in sorted(task.keys()):

  mytask = task[ t ] 
  cmd  = mytask[ 0 ]

  print("We have "+ t +"; "+ cmd)
  #r = call(["./scripts/submit_slurm.sh", t, cmd])
  #r=1

  #if r==0:
  #  i=i+1
  #  task[ t ] = [cmd, '', 0]
  #else:
  #  break


print "Prepared "+ str(len(task.keys())) +" jobs."

pickle.dump(task, open(filename, "wb") )
print "Saved in "+ filename

