#!/usr/bin/python
from subprocess import call
import pickle

filename = 'task_j0_13'

task={}
i=0
kk=0
r=0

#for k2 in range(-600,601,20):
#  for k1 in range(20,601,20):

for j2 in range(40,81,2):
  for j1 in range(40,101,2):
    #j2 = 30
    #k1=400
    sgn=1
  #for sgn in [0,1]:

#    t = "c3-"+str(k1)+'-'+str(k2)
    t = "jc1:"+str(j1)+":"+str(j2)

    #cmd = "main_u1kagome_scan23_c 12 1 0 1 0 1 1 100 "+ str(k) +" 50 "+ str(sgn) +" 0 0 200 300 8"
    #cmd = "main_u1kagome_scan23_c 12 1 0 1 0 1 1 100 50 50 "+ str(sgn) +" 0 0 "+ str(5*k) +" 300 10"
    #cmd = "main_u1compl 8 1 0 1 0 0 1 78 13 0 "+ str(k1) +" "+ str(k2) +" 0 10"
    #cmd = "main_u1compl 8 1 0 1 0 0 1 77 23 0 "+ str(k1) +" 0 0 10"
    #cmd = "main_u1kagome_scan12_c 12 1 0 1 1 0 1 100 "+ str(k) +" 50 "+ str(sgn) +" 0 300 0 0 8"
    #cmd = "main_he_kag 20 0 -20 "+ str(j2) +" "+ str(j3) +" 200"
    cmd = "main_he_kag 20 1 "+ str(j1) +" "+ str(j2) +" -40 200"
    #cmd = "main_he_kag 20 2 -20 "+ str(j2) +" "+ str(j3) +" 200"

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
