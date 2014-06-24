#!/usr/bin/env python

import pickle 

#create batch tasks that can be submitted to the grid later
filename = 'tasks_20.p'

tot=100
i=0
task={}

nn=0

#loop over the sign of xi3
#for x3sgn in [1, -1]:
  
  #loop over the sign of xi2
for x2sgn in [1, -1]:

    #loop over the phase (units of pi/600)
    for ph1 in range(0, 591, 50):

#      for xitot in range(2,99,4):
  
#        xi3 = int(x3sgn*(tot - xitot))
        xi3 = 0
        xitot = tot

#  cmd = "main_u1real_scan12 8 1 0 0 0 "+ str(xitot) +" -1 "+ str(xi3) +" 10",'',0
        #c-1, e2=1
        #cmd = "main_u1compl_scan12 8 1 0 0 1 2 1 "+ str(xitot) +" "+ str(x2sgn) +" "+ str(xi3) +" 0 "+ str(ph) +" 0 10"

        #c-2, e2=1 (second phase set to 0)
        #t_name = 'c2 '+str(x3sgn)+'_'+str(x2sgn)+'_'+str(ph)+'_'+str(xitot).zfill(2)
        #ph2=0
        #cmd = "main_u1compl_scan12 8 1 0 0 1 2 1 "+ str(xitot) +" "+ str(x2sgn) +" "+ str(xi3) +" "+ str(ph) +" 0 "+ str(ph2) +" 10"

        #c-3, e2=1
        for ph2 in range(0, 591, 50):
          if ph1==0 and ph2==0:
            continue
          if ph1 in [0, 150, 300, 450] and ph2 in [0, 150, 300, 450]:
            continue

          #t_name = 'c3:'+str(x3sgn)+'_'+str(x2sgn)+'_'+str(ph1)+'_'+str(ph2)+'_'+str(xitot).zfill(2)
          t_name = 'c3m:'+str(x2sgn)+'_'+str(ph1)+'_'+str(ph2)+'_'+str(xitot).zfill(2)
          #cmd = "main_u1compl_scan12 8 1 0 0 0 0 1 "+ str(xitot) +" "+ str(x2sgn) +" "+ str(xi3) +" "+ str(ph) +" "+ str(ph2) +" 0 10"
          cmd = "main_u1compl_scan12 8 1 0 1 0 0 1 "+ str(xitot) +" "+ str(x2sgn) +" "+ str(xi3) +" "+ str(ph1) +' '+ str(ph2) +" 0 10"

          #nn = nn+1
          #if nn>200:
          #  break
          task[ t_name ] = [cmd, '', 0, 'new'] #command, host, PID, status 
          i = i+1

pickle.dump(task, open(filename, "wb") )

for t in sorted(task.keys()):
  print t +': '+ task[t][0]

print str(len(task.keys()))
print "Ready to submit "+ str(i) +" jobs, stored in '"+ filename +"'"

