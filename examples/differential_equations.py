import sys
sys.path.append("..")

from sympy import Derivative, Symbol, Function, exp, Rational, log, \
    dsolve

import relativity


def eq1():
    r = Symbol("r")
    e = relativity.Rmn.dd(0,0)
    e = e.subs(relativity.nu(r), -relativity.lam(r))
    print dsolve(e, [relativity.lam(r)])

def eq2():
    r = Symbol("r")
    e = relativity.Rmn.dd(1,1)
    C = Symbol("CC")
    e = e.subs(relativity.nu(r), -relativity.lam(r))
    print dsolve(e, [relativity.lam(r)])

def eq3():
    r = Symbol("r")
    e = relativity.Rmn.dd(2,2)
    e = e.subs(relativity.nu(r), -relativity.lam(r))
    print dsolve(e, [relativity.lam(r)])

def eq4():
    r = Symbol("r")
    e = relativity.Rmn.dd(3,3)
    e = e.subs(relativity.nu(r), -relativity.lam(r))
    print dsolve(e, [relativity.lam(r)])

eq1()
eq2()
#eq3()
#eq4()
