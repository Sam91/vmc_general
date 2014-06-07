#!/usr/bin/env python

import pickle 

#create batch tasks that can be submitted to the grid later
filename = 'tasks_9.p'

tot=100
i=0
task={}

#loop over the sign of xi3
for x3sgn in [1, -1]:
  

  #loop over the sign of xi2
  for x2sgn in [1, -1]:

    #loop over the phase (units of pi/600)
    for ph in range(100, 501, 100):

      for xitot in range(2,99,2):
  
        xi3 = int(x3sgn*(tot - xitot))

#  task['task '+str(xitot)] = ["main_u1real_scan12 8 1 0 0 0 "+ str(xitot) +" -1 "+ str(xi3) +" 10",'',0]  
        task['j '+str(x3sgn)+'_'+str(x2sgn)+'_'+str(ph)+'_'+str(xitot)] = ["main_u1compl_scan12 8 1 0 0 1 2 1 "+ str(xitot) +" "+ str(x2sgn) +" "+ str(xi3) +" 0 "+ str(ph) +" 0 10",'',0]  

        i = i+1

pickle.dump(task, open(filename, "wb") )

for t in task.keys():
  print t +': '+ task[t][0]

print str(len(task.keys()))
print "Ready to submit "+ str(i) +" jobs, stored in '"+ filename +"'"

