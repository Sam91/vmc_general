#!/usr/bin/python

from subprocess import call
import pickle
import sys

task={}
i=0
kk=0
full = False

taskfile = str(sys.argv[1])
mesufile = str(sys.argv[2])

f = open(mesufile, "w")

#open tasklist
task = pickle.load( open(taskfile, "rb" ) )

for t in sorted(task.keys()):

  mytask = task[ t ]
  cmd  = mytask[0]
  srv  = mytask[1]
  stat = mytask[2]
  pid  = mytask[3]

  if stat==1 or stat == 'new':
    if not full:
      print task[ t ]
      f.write( '/home/bieri/vmc_general/bin/'+ cmd +' &\n' )
      i=i+1

      if i>=64:
        full=True
        break

      task[ t ] = [cmd, '', 'submitted', 0] #command, host, status, pid

  else:
    kk=kk+1

print "Submitted "+ str(i) +" new jobs; "+ str(len(task.keys()) - kk - i) +" still unsubmitted."

if i>0:
  pickle.dump(task, open(taskfile, "wb") )
  print "Saved in "+ taskfile

f.write( 'wait\n' )
f.close()

call( ['chmod','750',mesufile] )

