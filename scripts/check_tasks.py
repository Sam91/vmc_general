#!/usr/bin/env python

import sys
from batch_lib import *

taskfile = str(sys.argv[1])

tasklist = pickle.load( open(taskfile, "rb" ) )

checktasks( tasklist )

#submit(tasklist)

#pickle.dump(tasklist, open(taskfile, "wb"))

#killall(tasklist, 'cyclope')

#killjobs('sbieri', 'elektra')

#print tasklist


