#!/usr/bin/env python

import subprocess
import pickle
import time

import smtplib
from email.mime.text import MIMEText

servers = {'phenix','vega','cypher','christopher','maggott','fifi'}

all_servers = {'cyclope','thorgal','ghast','lockheed','christopher','mimic','dora','marrow','tapir','vega','beryl','cavallieri','wilbur','ember','colossus','yumi','pyro','bosons','trasses','callisto','sage','rox','cypher','iceberg','maggott','patch','mara','pongo','lisdal','sarabi','saravone','rosemayor','fifi','crush','cyborg','diego','choupette','bishop','bruce','solar','valerie','synch','droopy','grey','tornade','shadowcat','rovni','domino','duchesse','dazzler','clochette','biscotte','epervier','forge','magneto','khimaira','angel','peach','gambit','moonstar','havok','professorx','philos','phenix','prunelle','quinine','deb','polaris'}

#taking out haurele since Talbot is rebooting it
#taking out elektra for Frederic
#issue with 'akira' and 'donald'
#down: 'haurele', 'polaris'

servers2 = {'cyclope','thorgal','elektra',}

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
  r = 0.7

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
    nproc = int(sshexc(s, "ps -fe|grep -n 'bin/main_u1'|grep -v grep|wc -l"))

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
      print("Warn: No processes running on '"+ s +"'. You may check it")

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
def submit( tasklist ):

  tasks = sorted(tasklist.keys())

  #check if there is at least one job to submit
  no_new = True
  for t in tasks:
    if len(tasklist[ t ][0])>0:
      no_new = False
      break

  if no_new:
    print("No job to submit.")
    print "--------------------"

    checktasks( tasklist )
    return tasklist

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
    cmd = tasklist[ t ][0]
    srv = tasklist[ t ][1]

    if len(srv) > 0:
      #print "Task '"+ t +"' has already been submitted to '"+ srv +"'."
      #if srv=='akira' or srv=='donald':
      #  print("Exceptionally resubmitting "+ srv +" job")
      #else:
      continue
    waitingjob = waitingjob + 1

    [host, pid] = submit_job( cmd, nservers )

    if len(host)<1:
      break

    if pid==0:
      print("Warn: issue submitting to "+ host)
      continue

    tasklist[ t ][1] = host
    tasklist[ t ][2] = pid
    tasklist[ t ][3] = 'running'
    nsub = nsub+1
    newtask.append( t )

  #print tasklist
  print "--------------------"
  print "Submitted "+ str(nsub) +" jobs."
  print "--------------------"

  if nsub==0 and waitingjob==0:
    print('All pending jobs are submitted!')
    sendmail('All jobs submitted', 'All pending jobs for this tasklist have been submitted as of '+ dt)

  #save the updated servers list
  saveservers( nservers )

  if len(newtask)<1:
    newtask = None

  #check the running tasks
  checktasks( tasklist, newtask )

  return tasklist

#check the status of submitted tasks
def checktasks( tasklist, newtask = None ):

  done = 0
  running = 0
  notsub = 0
  unreachable = 0
  allkeys = sorted(tasklist.keys())

  tmsg = "Checking tasks on "+ time.strftime("%d/%m/%y %H:%M:%S")
  print( tmsg )

  for t in allkeys:
    srv = tasklist[t][1]
    if len(srv)<1:
      notsub = notsub + 1
      continue

    stat = tasklist[t][3]
    if stat=='done':
      done = done + 1
      continue
    if stat=='unreachable':
      unreachable = unreachable + 1
      continue
    if stat!='running':
      print("Warn: unknown state '"+ stat +"'")
      continue

    #check status for jobs marked as running
    pid = tasklist[t][2]
    r = isrunning(srv, pid)
    if r<0:
      print( "Server '"+ srv +"' unreachable.")
      unreachable = unreachable+1
      tasklist[t][3] = 'unreachable'
      continue
    if r==1:
      running = running + 1
    else:
      if r==0:
        done = done + 1
        tasklist[t][3] = 'done'
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

  return tasklist

#check if a job is running on a given server
def isrunning(srv, pid):
  if len(srv)<1 or pid<1:
    return -1

  out = sshexc(srv, "ps -fe|grep 'bin/main_u1'|grep ' "+ str(pid) +" '|grep -v grep|wc -l" )
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
      pid = tasklist[t][2]

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
          tasklist[t][2] == npid
          tasklist[t][3] == 'running'
          resubmitted = resubmitted + 1
        else:
          print("Invalid pid "+ str(npid) )

  print("Closed: "+ str(closed) +"\nResubmitted: "+ str(resubmitted))

#kill all submitted processed in the tasklist
def killall( tasklist ):

  for t in tasklist.keys():
    srv = tasklist[ t ][1]
    if len(srv)==0:
      continue
    pid = tasklist[ t ][2]

    print "Killing process '"+ str(pid) +"' on '"+ srv +"'."
    sshexc(srv, "kill -9 "+ str(pid) )

#kill the submitted tasks on a specific server
def killall( tasklist, server ):

  for t in tasklist.keys():
    srv = tasklist[ t ][1]
    if srv!=server:
      continue
    pid = tasklist[ t ][2]

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

