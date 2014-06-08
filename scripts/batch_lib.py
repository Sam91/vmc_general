#!/usr/bin/env python

import subprocess
import pickle
import time

import smtplib
from email.mime.text import MIMEText


servers = {'cyclope','thorgal','elektra','ghast','lockheed','christopher','mimic','dora','marrow','tapir','vega','beryl','cavallieri','wilbur','ember','colossus','yumi','pyro','bosons','trasses','callisto','sage','rox','cypher','iceberg','maggott','patch','mara','pongo','lisdal','sarabi','saravone','rosemayor','fifi','crush','cyborg','diego','choupette','bishop','bruce','solar','valerie','synch','droopy','grey','tornade','shadowcat','rovni','domino','duchesse','dazzler','clochette','biscotte','epervier','forge','magneto','khimaira','angel','peach','gambit','moonstar','havok','professorx','philos','phenix','prunelle','quinine',}

#down: 'haurele', 'polaris'

servers2 = {'cyclope','thorgal','elektra',}

def rexec( cmd ):
  pr = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
  output,stderr = pr.communicate()
  #print [output, stderr]
  return output

def sshexc(server, cmd):
  return rexec("ssh "+ server +" \""+ cmd +"\"" )


# Find the number of cpus and processors per server
def findservers():

  #fraction of CPUs to be used
  r = 0.8

  ntot = 0
  ctot = 0
  nused = 0

  print("Checking servers on "+ time.strftime("%d/%m/%y %H:%M:%S") )
  print("--------------------")

  nservers = {}
  for s in servers:

    #sshexc(s, "sleep 0") #dummy command to get the hostkey

    out = sshexc(s, "grep processor /proc/cpuinfo |wc -l")
    try:
      ncpu = int(out)
    except:
      print out
      continue
    nproc = int(sshexc(s, "ps -fe|grep -n 'bin/main_u1'|grep -v grep|wc -l"))

    #extract the idle percentage
    out = sshexc(s, "top -bn2 -d.2|grep 'Cpu(s)'|tail -1" )
    inx = out.find("id")
    try:
      idle = float( out[inx-6:inx] )
    except:
      print out
      continue

    print s + " has " + str(ncpu) + " CPUs and "+ str(nproc) + " running processes, idle = "+ str(idle)

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

#method to sumbit a list of task to all servers
def submit( tasklist ):

  #make sure we have an updated server list
  nservers = findservers()
  #nservers = loadservers()

  tasks = sorted(tasklist.keys())
  nsub = 0
  waitingjob = 0

  dt = time.strftime("%d/%m/%y %H:%M:%S")
  print("Submitting tasks on "+ dt )
  print("--------------------")

  for t in tasks:
    cmd = tasklist[ t ][0]
    srv = tasklist[ t ][1]

    if len(srv) > 0:
      #print "Task '"+ t +"' has already been submitted to '"+ srv +"'."
      continue
    waitingjob = waitingjob + 1

    host = gethost(nservers) #get host that has a free slot
    if( len(host)<1 ):
      break

    fcmd = 'vmc_general/scripts/local_submit.sh \''+ cmd +'\''

    out = sshexc(host, fcmd)
    try:
      pid = int(out)
    except:
      print("Unable to parse pid: "+ out)
      continue
    print host +"/"+ str(pid) +": "+ fcmd

    nservers[ host ][3] = nservers[ host ][3]-1
    #print nservers[ host ][3]

    tasklist[ t ][1] = host
    tasklist[ t ][2] = pid
    tasklist[ t ][3] = 'running'
    nsub = nsub+1

  #print tasklist
  print "--------------------"
  print "Submitted "+ str(nsub) +" jobs."
  print "--------------------"

  if nsub==0 and waitingjob==0:
    print('All pending jobs are submitted!')
    sendmail('All jobs submitted', 'All pending jobs for this tasklist have been submitted as of '+ dt)

  #save the updated servers list
  saveservers( nservers )

  #check the running tasks
  checktasks( tasklist )

  return tasklist

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

#check the status of submitted tasks
def checktasks( tasklist ):

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

    pid = tasklist[t][2]
    out = sshexc(srv, "ps -fe|grep 'bin/main_u1'|grep '"+ str(pid) +"'|grep -v grep|wc -l" )
    try:
      r = int(out)
    except:
      print( "Server '"+ srv +"' unreachable: "+ out )
      unreachable = unreachable+1
      tasklist[t][3] = 'unreachable'
      continue
    if r==1:
      running = running + 1
    else:
      if r==0:
        done = done + 1
        tasklist[t][3] = 'done'
      else:
        print "Warn: more than one jobs ("+ str(r) +") for "+ srv +"/"+ str(pid)

  n = len(allkeys)
  #msg = "Out of "+ str(n) +" tasks, "+ str(n-notsub) +" have been submitted; "+str(unreachable) +" are unreachable, "+ str(running) +" are running, and "+ str(done) +" are done."
  msg = "Total: "+ str(n) +"\nSubmitted: "+ str(n-notsub) +"\nUnreachable: "+ str(unreachable) +"\nRunning: "+ str(running) +"\nDone: "+ str(done)

  print( msg )
  sendmail('checktask report', tmsg +"\n"+ msg)

  #Send info mail when all jobs are done
  if done+unreachable>n:
    sendmail('All jobs done!', 'Happy :)' )

  return tasklist

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

