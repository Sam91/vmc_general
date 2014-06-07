#!/usr/bin/env python

import subprocess
import pickle

servers = {'cyclope','thorgal','elektra','ghast','lockheed','christopher','mimic','dora','marrow','tapir','vega','beryl','cavallieri','wilbur','ember','colossus','yumi','pyro','bosons','trasses','callisto','polaris','sage','rox','cypher','iceberg','maggott','patch','mara','pongo','lisdal','sarabi','saravone','rosemayor','fifi','crush','cyborg','diego','choupette','bishop','bruce','solar','valerie','synch','droopy','grey','tornade','shadowcat','rovni','domino','duchesse','dazzler','clochette','biscotte','epervier','forge','magneto','khimaira','angel','peach','gambit','moonstar','havok','professorx','philos','phenix','prunelle','quinine',}

#removing 'haurele' which is down.

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

  #fraction of idle CPUs to be used
  r = 0.8

  ntot = 0
  ctot = 0
  nused = 0

  nservers = {}
  for s in servers:

    sshexc(s, "sleep 0") #dummy command to get the hostkey

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
    ucpu = float(ncpu)*idle/100.
    if nproc>0:
      nsub =  max([0, int(round(r*ucpu)) - nproc])
    else:
      nsub = max([1,int(round(r*ucpu))])

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

  tasks = tasklist.keys()
  nsub = 0

  for t in tasks:
    cmd = tasklist[ t ][0]
    srv = tasklist[ t ][1]

    if len(srv) > 0:
      print "Task '"+ t +"' has already been submitted to '"+ srv +"'."
      continue

    host = gethost(nservers) #get host that has a free slot
    if( len(host)<1 ):
      break

    fcmd = 'vmc_general/scripts/local_submit.sh \''+ cmd +'\''

    pid = int(sshexc(host, fcmd))
    print host +"/"+ str(pid) +": "+ fcmd

    nservers[ host ][3] = nservers[ host ][3]-1
    #print nservers[ host ][3]

    tasklist[ t ][1] = host
    tasklist[ t ][2] = pid
    nsub = nsub+1

  #print tasklist
  print "Submitted "+ str(nsub) +" jobs."

  #save the updated servers list
  saveservers( nservers )
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
  allkeys = tasklist.keys()

  for t in allkeys:
    srv = tasklist[t][1]
    if len(srv)<1:
      notsub = notsub+1
      continue
    pid = tasklist[t][2]
    out = sshexc(srv, "ps -fe|grep 'bin/main_u1'|grep '"+ str(pid) +"'|grep -v grep|wc -l" )
    try:
      r = int(out)
    except:
      print( "Server '"+ srv +"' unreachable: "+ out )
      unreachable = unreachable+1
      continue
    if r==1:
      running = running+1
    else:
      if r==0:
        done = done +1
      else:
        print "Warn: more than one jobs ("+ str(r) +") for "+ srv +"/"+ str(pid)

  n = len(allkeys)
  print("Out of "+ str(n) +" tasks, "+ str(n-notsub) +" have been submitted; "+str(unreachable) +" are unreachable, "+ str(running) +" are running, and "+ str(done) +" are done.")


#Kill all jobs of a user on a specific server (be careful)
def killjobs(user, server):

  cmd = 'ps -u '+ user +' -o pid | grep -v "PID" | xargs kill -9'
  sshexc(server, cmd)

  print "Killed alll jobs of '"+ user +"' on '"+ server +"'"


