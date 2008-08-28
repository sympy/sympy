#!/usr/bin/python
#setup.py

import os

os.system('latex symbolicGA')
os.system('latex symbolicGA')
os.system('dvipdf symbolicGA')
os.system('rm *.aux *.out *.log *.dvi')
