
from utils import memoizer_immutable_args
from basic import Atom
from methods import ArithMeths, RelationalMeths

class NumberSymbol(ArithMeths, RelationalMeths, Atom):

    #@memoizer_immutable_args
    def __new__(cls):
        return object.__new__(cls)

class ImaginaryUnit(ArithMeths, RelationalMeths, Atom):

    #@memoizer_immutable_args
    def __new__(cls):
        return object.__new__(cls)

I = ImaginaryUnit()


class Exp1(NumberSymbol):

    pass

class Pi(NumberSymbol):
    def tostr(self):
        return "pi"

pi = Pi()

class GoldenRatio(NumberSymbol):

    pass

class EulerGamma(NumberSymbol):

    pass

class Catalan(NumberSymbol):

    pass

class NaN(NumberSymbol):

    pass

class Infinity(NumberSymbol):

    pass

class ComplexInfinity(NumberSymbol):

    pass
