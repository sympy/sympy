# This code will collect test info and store it in a database for viewing
#
# Collects:
#   Name of person who ran tests (looks for github name)
#   Python version used
#   Operating system used
#   Processor used
#   Run time for tests

import datetime
import os
import platform
import sqlite3
from sympy import *
import time
import timeit

# Collect user information and the current date
plat = platform.platform()
processor = platform.processor()
py_ver = platform.python_version()
date = datetime.date.today().isoformat()

# Attempt to find the users git name
process = os.popen('git config user.name')
username = process.read()
process.close()


# Create a function to determine average runtime
def average_time(command, runs=100000, setups=None):
    iterations = 30
    timeper = []
    for i in range(iterations):
        timeper.append(timeit.timeit(command, number=runs, setup=setups))

    return sum(timeper) / float(iterations)


# Establish a connection to the sqlite database
conn = sqlite3.connect('benchmarking.db')
c = conn.cursor()

# Lagrange test 1
# Simple Pendulum example
L1time = average_time('lagrange_1.lmethod()', 1000, setups='import lagrange_1')
values = (time.time(), username, date, plat, processor, py_ver, L1time)
c.execute("INSERT INTO lagrange_1 VALUES(%.6f, '%s', '%s', '%s', '%s', '%s',\
          %.6f)" % values)
conn.commit()

# Kane test 1
# Mass Spring Damper Docstring Example
K1time = average_time('kane_1.kmethod()', 1000, setups='import kane_1')
values = (time.time(), username, date, plat, processor, py_ver, K1time)
c.execute("INSERT INTO kane_1 VALUES(%.6f, '%s', '%s', '%s', '%s', '%s', %.6f)"
          % values)
conn.commit()

# Close the cursor and database connection properly
c.close()
conn.close()
