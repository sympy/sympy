#!/usr/bin/python

import sys
from IPython.core.display import display, Math

while True:
    line = sys.stdout.readline()
    #display(Math(line))
    sys.stderr.write('??? '+line.strip()+' ???')


