#!/usr/bin/python

import sys
from os import walk

f = []
for (dirpath, dirnames, filenames) in walk(mypath):
    f.extend(filenames)
    break

print str(filenames)
