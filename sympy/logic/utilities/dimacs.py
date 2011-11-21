"""For reading in DIMACS file format

www.cs.ubc.ca/~hoos/SATLIB/Benchmarks/SAT/satformat.ps

"""
from sympy.core import Symbol
from sympy.logic.boolalg import And, Or
import re

def load(s):
    clauses = []

    lines = s.split('\n')

    pComment = re.compile('c.*')
    pStats = re.compile('p\s*cnf\s*(\d*)\s*(\d*)')

    while len(lines) > 0:
        line = lines.pop(0)

        # Only deal with lines that aren't comments
        if not pComment.match(line):
            m = pStats.match(line)

            if not m:
                nums = line.rstrip('\n').split(' ')
                list = []
                for lit in nums:
                    if lit != '':
                        if int(lit) == 0: continue
                        num = abs(int(lit))
                        sign = True
                        if int(lit) < 0:
                            sign = False

                        if sign: list.append(Symbol("cnf_%s" % num))
                        else: list.append(~Symbol("cnf_%s" % num))

                if len(list) > 0:
                    clauses.append(Or(*list))

    return And(*clauses)

def load_file(location):
    s = open(location).read()
    return load(s)
