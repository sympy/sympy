import sys
sys.path.append("..")
from sympy.numerics import *
from sympy.numerics.utils_ import *
import math
from time import clock

def display_fraction(digits, skip=0, colwidth=10, columns=5):
    perline = colwidth * columns
    printed = 0
    for linecount in range((len(digits)-skip) // (colwidth * columns)):
        line = digits[skip+linecount*perline:skip+(linecount+1)*perline]
        for i in range(columns):
            print line[i*colwidth : (i+1)*colwidth],
        print ":", (linecount+1)*perline
        if (linecount+1) % 10 == 0:
            print
        printed += colwidth*columns
    rem = (len(digits)-skip) % (colwidth * columns)
    if rem:
        buf = digits[-rem:]
        s = ""
        for i in range(columns):
            s += buf[:colwidth].ljust(colwidth+1, " ")
            buf = buf[colwidth:]
        print s + ":", printed + colwidth*columns

def calculateit(func, base, n, tofile):
    Float.setprec(100)
    intpart = small_numeral(int(float(func())), base)
    if intpart == 0:
        skip = 0
    else:
        skip = len(intpart)
    Float.setprec(int(n*math.log(base,2))+10)
    print "Step 1 of 2: calculating binary value..."
    t = clock()
    a = func()
    step1_time = clock() - t
    print "Step 2 of 2: converting to specified base..."
    t = clock()
    d = bin_to_radix(a.man, -a.exp, base, n)
    d = fixed_to_str(d, base, n)
    step2_time = clock() - t
    print "\nWriting output...\n"
    if tofile:
        out_ = sys.stdout
        sys.stdout = tofile
    print "%i base-%i digits of pi:\n" % (n, base)
    print intpart, ".\n"
    display_fraction(d, skip, colwidth=10, columns=5)
    if tofile:
        sys.stdout = out_
    print "\nFinished in %f seconds (%f calc, %f convert)" % \
        ((step1_time + step2_time), step1_time, step2_time)

def interactive():
    print "Compute digits of pi with SymPy\n"
    base = input("Which base? (2-36, 10 for decimal) \n> ")
    digits = input("How many digits? (enter a big number, say, 10000)\n> ")
    tofile = raw_input("Output to file? (enter a filename, or just press enter\nto print directly to the screen) \n> ")
    if tofile:
        tofile = open(tofile, "w")
    global_options["verbose"] = True
    global_options["verbose_base"] = base
    calculateit(pi_float, base, digits, tofile)
    raw_input("\nPress enter to close this script.")

interactive()
