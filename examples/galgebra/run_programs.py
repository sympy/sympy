#!/usr/bin/python

import sys
from os import walk, system

f = []
for (dirpath, dirnames, filenames) in walk('./'):
    f.extend(filenames)
    break

for filename in f:
    if filename[-3:] == '.py':
        system('python ' + filename + '\n')
