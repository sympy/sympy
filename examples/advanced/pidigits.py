#!/usr/bin/env python

"""Pi digits example

Example shows arbitrary precision using mpmath with the
computation of the digits of pi.
"""

from mpmath import libmp, pi
from mpmath import functions as mpf_funs

import math
from time import clock
import sys


def display_fraction(digits, skip=0, colwidth=10, columns=5):
    """Pretty printer for first n digits of a fraction"""
    perline = colwidth * columns
    printed = 0
    for linecount in range((len(digits) - skip) // (colwidth * columns)):
        line = digits[skip + linecount*perline:skip + (linecount + 1)*perline]
        for i in range(columns):
            print(line[i*colwidth: (i + 1)*colwidth],)
        print(":", (linecount + 1)*perline)
        if (linecount + 1) % 10 == 0:
            print
        printed += colwidth*columns
    rem = (len(digits) - skip) % (colwidth * columns)
    if rem:
        buf = digits[-rem:]
        s = ""
        for i in range(columns):
            s += buf[:colwidth].ljust(colwidth + 1, " ")
            buf = buf[colwidth:]
        print(s + ":", printed + colwidth*columns)


def calculateit(func, base, n, tofile):
    """Writes first n base-digits of a mpmath function to file"""
    prec = 100
    intpart = libmp.numeral(3, base)
    if intpart == 0:
        skip = 0
    else:
        skip = len(intpart)
    print("Step 1 of 2: calculating binary value...")
    prec = int(n*math.log(base, 2)) + 10
    t = clock()
    a = func(prec)
    step1_time = clock() - t
    print("Step 2 of 2: converting to specified base...")
    t = clock()
    d = libmp.bin_to_radix(a.man, -a.exp, base, n)
    d = libmp.numeral(d, base, n)
    step2_time = clock() - t
    print("\nWriting output...\n")
    if tofile:
        out_ = sys.stdout
        sys.stdout = tofile
    print("%i base-%i digits of pi:\n" % (n, base))
    print(intpart, ".\n")
    display_fraction(d, skip, colwidth=10, columns=5)
    if tofile:
        sys.stdout = out_
    print("\nFinished in %f seconds (%f calc, %f convert)" % \
        ((step1_time + step2_time), step1_time, step2_time))


def interactive():
    """Simple function to interact with user"""
    print("Compute digits of pi with SymPy\n")
    base = input("Which base? (2-36, 10 for decimal) \n> ")
    digits = input("How many digits? (enter a big number, say, 10000)\n> ")
    tofile = raw_input("Output to file? (enter a filename, or just press enter\nto print directly to the screen) \n> ")
    if tofile:
        tofile = open(tofile, "w")
    calculateit(pi, base, digits, tofile)


def main():
    """A non-interactive runner"""
    base = 16
    digits = 500
    tofile = None
    calculateit(pi, base, digits, tofile)

if __name__ == "__main__":
    interactive()
