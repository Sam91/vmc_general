#!/usr/bin/env python

import subprocess
import socket
import nmap
import paramiko

from batch_lib import *

ssh = paramiko.SSHClient()
ssh.load_host_keys("/users/invites/sbieri/.ssh/my_known_hosts")
ssh.set_missing_host_key_policy( paramiko.AutoAddPolicy() )

nm = nmap.PortScanner()
nm.scan(hosts='134.157.8.0/24', arguments='-n -sP')

host_list=[(x,nm[x]['status']['state']) for x in nm.all_hosts()]

hosts={}

for a in host_list:

  ip = str(a[0])

  try:
    [addr,s,s] = socket.gethostbyaddr( ip )
  except:
    addr = 'unknown'
  hosts[ip] = addr

  try:
    out = ssh.connect( ip )
  except:
    out = 'error'
    continue

  host = addr[0:addr.find('.')]

  found = False
  for s in servers:
    if s==host:
      found=True
      break

  if not found:
    print [ip, addr, out]


