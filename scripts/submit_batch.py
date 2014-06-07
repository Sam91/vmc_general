#!/usr/bin/env python

from batch_lib import *

taskfile = 'tasks_tst.p'

tasklist = pickle.load( open(taskfile, "rb" ) )

submit(tasklist)

pickle.dump(tasklist, open(taskfile, "wb"))

#checktasks( tasklist )

#killall(tasklist, 'cyclope')

#killjobs('sbieri', 'elektra')

#print tasklist


