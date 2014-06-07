#!/usr/bin/env python

import sys
from batch_lib import *

#taskfile = 'tasks_tst.p'
taskfile = str(sys.argv[1])

#open tasklist
tasklist = pickle.load( open(taskfile, "rb" ) )

submit(tasklist)

#save updated tasklist
pickle.dump(tasklist, open(taskfile, "wb"))

#checktasks( tasklist )

#killall(tasklist, 'cyclope')

#killjobs('sbieri', 'elektra')

#print tasklist


