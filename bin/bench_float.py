# basic benchmark for the Float class

from time import clock
from random import *
from sympy.numerics import *
from sympy.numerics.functions import *
from sympy.numerics.functions2 import *

w = """
def f(N):
    xs = []
    seed(1234)
    for i in range(N):
        x = Float(random() * 2.0**randint(-10, 10)) ** 0.5
        y = Float(random() * 2.0**randint(-10, 10)) ** 0.5
        x *= choice([1, -1])
        y *= choice([1, -1])
        xs.append((x, y))
    t1 = clock()
    for x, y in xs:
        OP
    t2 = clock()
    return int(N/(t2-t1))

tests.append(("OP", f))
"""

tests = []
def deftest(op):
    exec w.replace("OP", op)

atests = ["x + y", "x - y", "x * y", "x / y", "x == y",
  "x < y", "abs(x)", "abs(x)**0.5", "exp(x)", "sin(x)",
  "tan(x)", "atan(x)"]

for test in atests:
    deftest(test)

precs = [15, 30, 100, 500, 1000]

def runtests():
    print "\n prec (dps) =",
    for prec in precs:
        print "%7s" % prec,
    print
    print "-" * 75
    for name, test in tests:
        print ("%12s" % name), ":",
        for prec in precs:
            Float.setdps(prec)
            print "%7s" % test(max(50, 5000//prec)),
            Float.revert()
        print

print "\nFloat timings (operations / second)"

runtests()
