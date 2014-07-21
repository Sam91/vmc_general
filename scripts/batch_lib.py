#!/usr/bin/env python

import subprocess
import pickle
import time

import smtplib
from email.mime.text import MIMEText

#servers = {'mimic','iceberg','gambit','tornade','cypher','bishop','rovni','trasses','droopy','yumi','domino','philos', 'deb', 'ghast', 'phenix','lockheed','maggott','fifi','cavallieri','marrow'}

#servers = {'ghast','lockheed','mimic','marrow','tapir','vega','beryl','cavallieri','ember','colossus','yumi','pyro','bosons','trasses','callisto','sage','rox','cypher','iceberg','maggott','patch','mara','pongo','lisdal','sarabi','saravone','rosemayor','fifi','crush','cyborg','diego','choupette','bishop','bruce','solar','valerie','synch','droopy','grey','tornade','shadowcat','rovni','domino','duchesse','dazzler','clochette','biscotte','epervier','forge','magneto','khimaira','angel','peach','gambit','moonstar','havok','professorx','philos','phenix','prunelle','quinine','deb','polaris'}

servers = ['phenix','philos','maggott','iceberg','fifi','bishop','droopy','vega','clochette','choupette','christopher','bosons','polaris','magneto','professorx','solar','crush','trasses','forge','cypher']

#taking out 'christopher' because it is slow
#taking out 'haurele' since Talbot is rebooting it
#taking out 'elektra' for Frederic
#issue with 'akira' and 'donald'
#mimeto keeps rebooting
#down: 'haurele', 'polaris'

f256 = {'cyclope','dora','thorgal','wilbur'}

servers2 = {'elektra'}

#execute a command in the shell
def rexec( cmd ):
  #print( cmd )
  pr = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
  output,stderr = pr.communicate()
  #print [output, stderr]
  return output

#execute a remote command vie ssh
def sshexc(server, cmd):
  return rexec("ssh -o BatchMode=yes "+ server +" \""+ cmd +"\"" )

# Find the number of cpus and processes per server
def findservers():

  #fraction of CPUs to be used
  r = 0.8

  ntot = 0
  ctot = 0
  nused = 0

  print("Parsing servers on "+ time.strftime("%d/%m/%y %H:%M:%S") )
  print("--------------------")

  nservers = {}
  for s in sorted(servers):

    #sshexc(s, "sleep 0") #dummy command to get the hostkey

    out = sshexc(s, "grep processor /proc/cpuinfo |wc -l")
    try:
      ncpu = int(out)
    except:
      print out
      continue
    nproc = int(sshexc(s, "ps -fe|grep -n 'bin/main_'|grep -v grep|grep -v srun|wc -l"))

    #extract the idle percentage
    out = sshexc(s, "top -bn2 -d.2|grep 'Cpu(s)'|tail -1" ).replace('%',' ')
    inx = out.find("id")
    try:
      idle = float( out[inx-6:inx] )
    except:
      print out
      continue

    print s + " has " + str(ncpu) + " CPUs and "+ str(nproc) + " running processes, idle = "+ str(idle)

    if nproc==0 and (.8*idle/100.>1.-r):
      print("==> Warn: No processes running on '"+ s +"'. You may check it")

    #calculate the number of jobs we want to sumbit to this server
    #theoretical maximum
    nmax = max([0, int(round(r*float(ncpu)))-nproc])
    if nmax>0:
      nsub = max([0, int(round((idle/100.+r-1.)*float(ncpu)))])
    else:
      nsub = 0
    if nsub > nmax:
      nsub = nmax

    nservers[ s ] = [ncpu, nproc, idle, nsub]

    ntot = ntot + nsub
    ctot = ctot + ncpu
    nused = nused + nproc

  print "------------"
  print "We can potentially run "+ str(ntot) +" more jobs; "+ str(nused) +" currently running. CPUs: "+ str(ctot)
  print "------------"

  saveservers( nservers )
  return nservers

#save the servers list to file
def saveservers(s):
  pickle.dump(s, open("nservers.p", "wb") )

#load the server list from a file (be careful: it may be outdated)
def loadservers():
  return pickle.load( open("nservers.p", "rb" ) )

#return a host that has at least one free slot
def gethost(nservers):
  for s in nservers.keys():
    if nservers[ s ][3] > 0:
      return s
  print "No free host found..."
  return ''

#submit a job via slurm
def submit_slurmjob( job, cmd ):

  fcmd = 'scripts/submit_slurm.sh '+ job +" '"+ cmd +"'"
  out = rexec( fcmd )

  return out

def submit_slurm( tasklist ):

  tasks = sorted( tasklist.keys() )

  nsub = 0
  waitingjob = 0

  dt = time.strftime("%d/%m/%y %H:%M:%S")
  print("Submitting tasks on "+ dt )
  print("--------------------")

  newtask=[]
  for t in tasks:
    cmd = tasklist[ t ][0]
    srv = tasklist[ t ][1]
    sts = tasklist[ t ][2]

    if sts == 0:
      continue
    waitingjob = waitingjob + 1

    out = submit_slurmjob( t, cmd )

    #if len(host)<1:
    #  break

#    if pid==0:
#      print("Warn: issue submitting to "+ host)
#      continue
    print out


    tasklist[ t ][1] = out
    #tasklist[ t ][2] = pid
    #tasklist[ t ][3] = 'running'
    tasklist[ t ][2] = 'submitted'
    nsub = nsub+1
    newtask.append( t )

  #print tasklist
  print "--------------------"
  print "Submitted "+ str(nsub) +" jobs."
  print "--------------------"

#  if nsub==0 and waitingjob==0:
#    print('All pending jobs are submitted!')
#    sendmail('All jobs submitted', 'All pending jobs for this tasklist have been submitted as of '+ dt)

  #save the updated servers list
  #saveservers( nservers )

  if len(newtask)<1:
    newtask = None

  #check the running tasks
#  checktasks( tasklist, newtask )

  return tasklist


#submit a job to a server
def submit_job( cmd, nservers ):

    host = gethost(nservers) #get host that has a free slot
    if( len(host)<1 ):
      return ['', 0]

    fcmd = 'vmc_general/scripts/local_submit.sh \''+ cmd +'\''

    out = sshexc(host, fcmd)
    try:
      pid = int(out)
    except:
      print("Unable to parse pid: "+ out)
      return [host, 0]
    print host +"/"+ str(pid) +": "+ fcmd

    nservers[ host ][3] = nservers[ host ][3]-1

    return [host, pid]

#method to sumbit a list of task to all servers
def submit( task ):

  tasks = sorted(task.keys())

  #check if there is at least one job to submit
  no_new = True
  for t in tasks:
    if task[t][2]=='new':
      no_new = False
      break

  if no_new:
    print("No job to submit.")
    print "--------------------"

    checktasks( task )
    return task

  #make sure we have an updated server list
  nservers = findservers()
  #nservers = loadservers()

  nsub = 0
  waitingjob = 0

  dt = time.strftime("%d/%m/%y %H:%M:%S")
  print("Submitting tasks on "+ dt )
  print("--------------------")

  newtask=[]
  for t in tasks:

    cmd = task[ t ][0] #command
    srv = task[ t ][1] #server
    sts = task[ t ][2] #status

    if sts=='done' or sts=='submitted' or sts=='running':
      #print "Task '"+ t +"' has already been submitted to '"+ srv +"'."
      #if srv=='akira' or srv=='donald':
      #if srv=='mimeto':
      #  print("Exceptionally resubmitting "+ srv +" job")
      #else:
        continue
    waitingjob = waitingjob + 1

    if sts=='unreachable':
    #if sts=='unreachable' and len(srv)<1:
      continue

    #print "Submitting job in state '"+ str(sts) +"'"
    print task[t]

    [host, pid] = submit_job( cmd, nservers )

    if len(host)<1:
      break

    if pid==0:
      print("Warn: issue submitting to "+ host)
      continue

    task[ t ][1] = host
    task[ t ][2] = 'running'
    task[ t ][3] = pid
    nsub = nsub+1
    newtask.append( t )

  #print tasklist
  print "--------------------"
  print "Submitted "+ str(nsub) +" jobs."
  print "--------------------"

  if nsub==0 and waitingjob==0:
    print('All pending jobs are submitted!')
    #sendmail('All jobs submitted', 'All pending jobs for this tasklist have been submitted as of '+ dt)

  #save the updated servers list
  saveservers( nservers )

  if len(newtask)<1:
    newtask = None

  #check the running tasks
  checktasks( task, newtask )

  return task

#check the status of submitted tasks
def checktasks( task, newtask = None ):

  done = 0
  running = 0
  notsub = 0
  unreachable = 0
  allkeys = sorted(task.keys())

  tmsg = "Checking tasks on "+ time.strftime("%d/%m/%y %H:%M:%S")
  print( tmsg )

  for t in allkeys:
    srv = task[t][1]
    if len(srv)<1:
      notsub = notsub + 1
      continue

    stat = task[t][2]
    if stat=='done' or stat=='submitted' or stat==0:
      done = done + 1
      continue
    if stat=='unreachable':
      unreachable = unreachable + 1
      continue
    if len(srv)<1 or stat=='new':
      continue

    #check status for jobs marked as running
    if stat!='running':
      print("Warn: unknown state '"+ stat +"'")
      continue

    pid = task[t][3]
    r = isrunning(srv, pid)
    if r<0:
      print( "Server '"+ srv +"' unreachable.")
      unreachable = unreachable+1
      task[t][2] = 'unreachable'
      continue
    if r==1:
      running = running + 1
    else:
      if r==0:
        done = done + 1
        task[t][2] = 'done'
        if newtask is not None:
          if t in newtask:
            print("Warn: New job was NOT found: "+ srv +"/"+ str(pid) )
      else:
        print "Warn: more than one job ("+ str(r) +") for "+ srv +"/"+ str(pid)

  if newtask is None:
    nnew = 0
  else:
    nnew = len(newtask)

  n = len(allkeys)
  #msg = "Out of "+ str(n) +" tasks, "+ str(n-notsub) +" have been submitted; "+str(unreachable) +" are unreachable, "+ str(running) +" are running, and "+ str(done) +" are done."
  msg = "Total: "+ str(n) +"\nSubmitted: "+ str(n-notsub) +"\nNew: "+ str(nnew) +"\nUnreachable: "+ str(unreachable) +"\nRunning: "+ str(running) +"\nDone: "+ str(done)

  print( msg )
  sendmail('checktask report', tmsg +"\n"+ msg)

  #Send info mail when all jobs are done
  if done+unreachable>n:
    sendmail('All jobs done!', 'Happy :)' )

  return task

#check if a job is running on a given server
def isrunning(srv, pid):
  if len(srv)<1 or pid<1:
    return -1

  out = sshexc(srv, "ps -fe|grep 'bin/main_'|grep ' "+ str(pid) +" '|grep -v grep|wc -l" )
  r = 0
  try:
    r = int(out)
  except:
    r = -1

  if r>1:
    print("Warn: more than one processes; "+ srv +"/"+ str(pid))

  return r

#go through unreachable jobs and try to fix them
def fixunreachable( tasklist ):

  tasks = sorted(tasklist.keys())

  no_unreachable = True
  for t in tasks:
    if tasklist[t][3] == 'unreachable':
      no_unreachable = False
      break

  if no_unreachable:
    return

  #get an updated server list
  nservers = findservers()
  #nservers = loadservers()

  closed=0
  resubmitted=0
  for t in tasks:
    if tasklist[t][3] == 'unreachable':

      srv = tasklist[t][1]
      pid = tasklist[t][3]

      r = isrunning(srv, pid)

      reachable = True
      if r<0:
        print( "Server '"+ srv +"' still seems unreachable; "+ out )
        reachable = False

      if reachable:
        print("Server '"+ srv +"' seems to be reachable again.")
        if r==0:
          print("Setting task to done.")
          tasklist[t][3] = 'done'
          closed = closed + 1
        continue

      else:
        print('Resubmitting job...')
        [host, npid] = submit_job( tasklist[t][0], nservers )
        if len(host)<1:
          break
        if npid>0:
          tasklist[t][1] == host
          tasklist[t][2] == 'running'
          tasklist[t][3] == npid
          resubmitted = resubmitted + 1
        else:
          print("Invalid pid "+ str(npid) )

  print("Closed: "+ str(closed) +"\nResubmitted: "+ str(resubmitted))

#kill all submitted processed in the tasklist
def kill( tasklist ):

  for t in sorted(tasklist.keys()):
    srv = tasklist[ t ][1]
    if len(srv)==0:
      continue
    pid = tasklist[ t ][3]

    print "Killing process '"+ str(pid) +"' on '"+ srv +"'."
    sshexc(srv, "kill -9 "+ str(pid) )

#kill the submitted tasks on a specific server
def killall( tasklist, server ):

  for t in tasklist.keys():
    srv = tasklist[ t ][1]
    if srv!=server:
      continue
    pid = tasklist[ t ][3]

    print "Killing process '"+ str(pid) +"' on '"+ srv +"'."
    sshexc(srv, "kill -9 "+ str(pid) )

#Kill all jobs of a user on a specific server (be careful)
def killjobs(user, server):

  cmd = 'ps -u '+ user +' -o pid | grep -v "PID" | xargs kill -9'
  sshexc(server, cmd)

  print "Killed alll jobs of '"+ user +"' on '"+ server +"'"

#send a mail to predefined user
def sendmail(subj, txt):

  msg = MIMEText( txt )

  me = 'sbieri@lptl.jussieu.fr'
  you = 'sbieri@lptl.jussieu.fr'

  msg['Subject'] = subj
  msg['From'] = me
  msg['To'] = you

  s = smtplib.SMTP('bambi.lptl.jussieu.fr')
  s.sendmail(me, [you], msg.as_string())
  s.quit()

