#!/usr/bin/env python

#takes a task file as argument and updates the status of each job

import sys
from batch_lib import *

taskfile = str(sys.argv[1])

tasklist = pickle.load( open(taskfile, "rb" ) )

fixunreachable( tasklist )

pickle.dump(tasklist, open(taskfile, "wb"))

