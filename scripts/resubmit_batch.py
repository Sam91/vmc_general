#!/usr/bin/python

from subprocess import call
import pickle
import sys

task={}
i=0
kk=0
full = False

taskfile = str(sys.argv[1])
batchfile = str(sys.argv[2])

#open tasklist
task = pickle.load( open(taskfile, "rb" ) )

f = open(batchfile, "a" )

for t in sorted(task.keys()):

  mytask = task[ t ]
  cmd  = mytask[0]
  srv  = mytask[1]
  pid  = mytask[2]
  stat = mytask[3]

  if stat==1 or stat == 'new':
    if not full:
      print task[ t ]

      i=i+1
      if i>64:
        full=True
        continue

      cmd = "/users/bieri/vmc_general/bin/"+ cmd +" &\n"
      f.write( cmd )

      task[ t ] = [cmd, '', 0, 'submitted'] #command, host, status, pid
    else:
      print task[ t ]
    #  task[ t ] = [cmd, srv, stat, pid]
  else:
    if stat == 'new':
      print tark[ t ]
    #task[ t ] = [cmd, srv, stat, pid]
    kk=kk+1

print "Submitted "+ str(i) +" new jobs; "+ str(len(task.keys()) - kk - i) +" still unsubmitted."

if i>0:
  pickle.dump(task, open(taskfile, "wb") )
  print "Saved in "+ taskfile

f.close()

