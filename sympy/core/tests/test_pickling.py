import pickle
from sympy.utilities.pytest import XFAIL

from sympy.core.assumptions import AssumeMeths
from sympy.core.basic import Atom, Basic, BasicMeta, BasicType,\
                             ClassesRegistry, Singleton, SingletonFactory
from sympy.core.symbol import Dummy, Symbol, Temporary, Wild
from sympy.core.numbers import Catalan, ComplexInfinity, EulerGamma, Exp1,\
                               GoldenRatio, Half, ImaginaryUnit, Infinity,\
                               Integer, NaN, NegativeInfinity,  NegativeOne,\
                               Number, NumberSymbol, One, Pi, Rational, Real,\
                               Zero
from sympy.core.relational import Equality, Inequality, Relational,\
                                  StrictInequality, Unequality
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.core.methods import ArithMeths, NoArithMeths, NoRelMeths, RelMeths

from sympy.core.function import Derivative, Function, FunctionClass, Lambda,\
                                WildFunction
from sympy.core.interval import Interval
from sympy.core.multidimensional import vectorize
from sympy.core.cache import Memoizer
from sympy.core.ast_parser import SymPyParser, SymPyTransformer


def check(a):
    b = pickle.loads(pickle.dumps(a))
    d1 = dir(a)
    d2 = dir(b)
    assert d1==d2

    def c(a,b,d):
        for i in d:
            if not hasattr(a,i):
                continue
            attr = getattr(a,i)
            if not hasattr(attr, "__call__"):
                assert hasattr(b,i), i
                assert getattr(b,i)==attr
    c(a,b,d1)
    c(b,a,d2)

def test_core_assumptions():
    for c in (AssumeMeths, AssumeMeths()):
        check(c)

def test_core_basic():
    for c in (Atom, Atom(), Basic, Basic(), BasicMeta, BasicMeta("test"),
              BasicType, BasicType("test"), ClassesRegistry, ClassesRegistry(),
              Singleton, Singleton(), SingletonFactory, SingletonFactory()):
        check(c)

def test_core_symbol():
    for c in (Dummy, Dummy("x", False), Symbol, Symbol("x", False),
              Temporary, Temporary(), Wild, Wild("x")):
        check(c)

def test_core_numbers():
    for c in (Catalan, Catalan(), ComplexInfinity, ComplexInfinity(),
              EulerGamma, EulerGamma(), Exp1, Exp1(), GoldenRatio, GoldenRatio(),
              Half, Half(), ImaginaryUnit, ImaginaryUnit(), Infinity, Infinity(),
              Integer, Integer(2), NaN, NaN(), NegativeInfinity,
              NegativeInfinity(), NegativeOne, NegativeOne(), Number, Number(15),
              NumberSymbol, NumberSymbol(), One, One(), Pi, Pi(), Rational,
              Rational(1,2), Real, Real("1.2"), Zero, Zero()):
        check(c)

def test_core_relational():
    x = Symbol("x")
    y = Symbol("y")
    for c in (Equality, Equality(x,y), Inequality, Inequality(x,y), Relational,
              Relational(x,y), StrictInequality, StrictInequality(x,y), Unequality,
              Unequality(x,y)):
        check(c)

def test_core_add():
    x = Symbol("x")
    for c in (Add, Add(x,4)):
        check(c)

def test_core_mul():
    x = Symbol("x")
    for c in (Mul, Mul(x,4)):
        check(c)

def test_core_power():
    x = Symbol("x")
    for c in (Pow, Pow(x,4)):
        check(c)

def test_core_methods():
    for c in (ArithMeths, ArithMeths(), NoArithMeths, NoArithMeths(),
              NoRelMeths, NoRelMeths(), RelMeths, RelMeths()):
        check(c)

def test_core_function():
    x = Symbol("x")
    for f in (Derivative, Derivative(x), Function, FunctionClass, Lambda,\
              WildFunction):
        check(f)

@XFAIL
def test_core_dynamicfunctions():
    # This fails because f is assumed to be a class at sympy.basic.function.f
    f = Function("f")
    check(f)

def test_core_interval():
    for c in (Interval, Interval(0,2)):
        check(c)

def test_core_multidimensional():
    for c in (vectorize, vectorize(0)):
        check(c)

@XFAIL
def test_core_cache():
    for c in (Memoizer, Memoizer()):
        check(c)

@XFAIL
def test_core_astparser():
    # This probably fails because of importing the global sympy scope.
    for c in (SymPyParser, SymPyParser(), SymPyTransformer,
              SymPyTransformer({},{})):
        check(c)
