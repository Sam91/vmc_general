#!/usr/bin/env python

import pickle 

#create batch tasks that can be submitted to the grid later
filename = 'tasks_tst.p'

tot=100
i=0
task={}

for xitot in range(2,98,40):

  xi3 = -(tot - xitot)

  task['task '+str(xitot)] = ["main_u1real_scan12 8 1 0 0 0 "+ str(xitot) +" -1 "+ str(xi3) +" 10",'',0]  

  i = i+1

pickle.dump(task, open(filename, "wb") )

for t in task.keys():
  print t +': '+ task[t][0]

print "Ready to submit "+ str(i) +" jobs, stored in '"+ filename +"'"

