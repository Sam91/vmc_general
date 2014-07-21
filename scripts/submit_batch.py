#!/usr/bin/env python

import sys
from batch_lib import *

#taskfile = '/users/invites/sbieri/vmc_general/tasks_20.p'
taskfile = str(sys.argv[1])

#open tasklist
task = pickle.load( open(taskfile, "rb" ) )

submit( task )
#kill( task )

#save updated tasklist
pickle.dump(task, open(taskfile, "wb"))

#submit a second taskfile
#taskfile = '/users/invites/sbieri/vmc_general/tasks_12.p'
#tasklist = pickle.load( open(taskfile, "rb" ) )
#submit(tasklist)
#pickle.dump(tasklist, open(taskfile, "wb"))


#killall(tasklist, 'cyclope')

#killjobs('sbieri', 'elektra')

#print tasklist

