#!/usr/bin/python

from subprocess import call
import pickle
import sys

task={}
i=0
kk=0
full = False

taskfile = str(sys.argv[1])

#open tasklist
task = pickle.load( open(taskfile, "rb" ) )

for t in sorted(task.keys()):

  mytask = task[ t ]
  cmd = mytask[0]
  srv = mytask[1]
  stat = mytask[2]
  pid = mytask[3]
  #pid = 0

  if stat==1 or stat == 'new':
    if not full:
      print task[ t ]
      r = call(["/users/invites/sbieri/vmc_general/scripts/submit_slurm.sh", t, cmd])

      if r==0:
        i=i+1
      else:
        full=True
        continue

      task[ t ] = [cmd, '', 'submitted', 0] #command, host, status, pid
    #else:
    #  task[ t ] = [cmd, srv, stat, pid]
  else:
    #task[ t ] = [cmd, srv, stat, pid]
    kk=kk+1

print "Submitted "+ str(i) +" new jobs; "+ str(len(task.keys()) - kk - i) +" still unsubmitted."

if i>0:
  pickle.dump(task, open(taskfile, "wb") )
  print "Saved in "+ taskfile

