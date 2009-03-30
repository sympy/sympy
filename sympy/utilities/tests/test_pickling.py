import copy
import pickle
import types
from sympy.utilities.pytest import XFAIL

from sympy.core.assumptions import AssumeMeths
from sympy.core.basic import Atom, Basic, BasicMeta, BasicType,\
        ClassesRegistry, SingletonFactory
from sympy.core.symbol import Dummy, Symbol, Temporary, Wild
from sympy.core.numbers import Catalan, ComplexInfinity, EulerGamma, Exp1,\
        GoldenRatio, Half, ImaginaryUnit, Infinity, Integer, NaN,\
        NegativeInfinity,  NegativeOne, Number, NumberSymbol, One, Pi,\
        Rational, Real, Zero
from sympy.core.relational import Equality, Inequality, Relational,\
        StrictInequality, Unequality
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.core.function import Derivative, Function, FunctionClass, Lambda,\
        WildFunction
from sympy.core.interval import Interval
from sympy.core.multidimensional import vectorize
from sympy.core.cache import Memoizer
#from sympy.core.ast_parser import SymPyParser, SymPyTransformer

from sympy import symbols


def check(a):
    """ Check that pickling and copying round-trips.
    """
    for protocol in [0, 1, 2, copy.copy, copy.deepcopy]:
        if callable(protocol):
            if isinstance(a, BasicType):
                # Classes can't be copied, but that's okay.
                return
            b = protocol(a)
        else:
            b = pickle.loads(pickle.dumps(a, protocol))

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


#================== core =========================

def test_core_assumptions():
    for c in (AssumeMeths, AssumeMeths()):
        check(c)

def test_core_basic():
    for c in (Atom, Atom(), Basic, Basic(), BasicMeta, BasicMeta("test"),
              BasicType, BasicType("test"), ClassesRegistry, ClassesRegistry(),
              SingletonFactory, SingletonFactory()):
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

# This doesn't have to be pickable.
#@XFAIL
#def test_core_astparser():
#    # This probably fails because of importing the global sympy scope.
#    for c in (SymPyParser, SymPyParser(), SymPyTransformer,
#              SymPyTransformer({},{})):
#        check(c)


#================== functions ===================
from sympy.functions import (Piecewise, lowergamma, acosh,
        chebyshevu, chebyshevt, ln, chebyshevt_root, Binomial, legendre,
        Heaviside, Dij, factorial, bernoulli, coth, tanh, assoc_legendre, sign,
        arg, asin, DiracDelta, re, rf, abs, uppergamma, binomial, sinh, Ylm,
        oo, cos, cot, acos, acot, gamma, bell, hermite, harmonic,
        LambertW, zeta, log, Factorial, pi, asinh, acoth, Zlm,
        cosh, dirichlet_eta, Eijk, loggamma, erf, max_, ceiling, im, fibonacci,
        conjugate, tan, chebyshevu_root, floor, atanh, nan, sqrt, zoo, min_,
        RisingFactorial, sin, E, atan, I, ff, FallingFactorial, lucas, atan2,
        polygamma, exp)

def test_functions():
    zero_var = (pi, oo, nan, zoo, E, I)
    one_var = (acosh, ln, Heaviside, Dij, factorial, bernoulli, coth, tanh,
            sign, arg, asin, DiracDelta, re, abs, sinh, cos, cot, acos, acot,
            gamma, bell, harmonic, LambertW, zeta, log, Factorial, asinh,
            acoth, cosh, dirichlet_eta, loggamma, erf, ceiling, im, fibonacci,
            conjugate, tan, floor, atanh, sin, atan, lucas, exp)
    two_var = (rf, ff, lowergamma, chebyshevu, chebyshevt, binomial, max_,
            min_, atan2, polygamma, hermite, legendre, uppergamma)
    x, y, z = symbols("x y z")
    others = (chebyshevt_root, chebyshevu_root, Eijk(x, y, z),
            Piecewise( (0, x<-1), (x**2, x<=1), (x**3, True)),
            assoc_legendre)
    for a in zero_var:
        check(a)
    for cls in one_var:
        check(cls)
        c = cls(x)
        check(c)
    for cls in two_var:
        check(cls)
        c = cls(x, y)
        check(c)
    for cls in others:
        check(cls)

#================== geometry ====================
from sympy.geometry.entity import GeometryEntity
from sympy.geometry.point import Point
from sympy.geometry.ellipse import Circle, Ellipse
from sympy.geometry.line import Line, LinearEntity, Ray, Segment
from sympy.geometry.polygon import Polygon, RegularPolygon, Triangle

def test_geometry():
    p1 = Point(1,2)
    p2 = Point(2,3)
    p3 = Point(0,0)
    p4 = Point(0,1)
    for c in (GeometryEntity, GeometryEntity(), Point, p1, Circle, Circle(p1,2),
              Ellipse, Ellipse(p1,3,4), Line, Line(p1,p2), LinearEntity,
              LinearEntity(p1,p2), Ray, Ray(p1,p2), Segment, Segment(p1,p2),
              Polygon, Polygon(p1,p2,p3,p4), RegularPolygon, RegularPolygon(p1,4,5),
              Triangle):
        # XXX: Instance of Triangle hangs because of hasattr in check().
        # Triangle(p1,p2,p3)
        check(c)
        pass

#================== integrals ====================
from sympy.integrals.integrals import Integral

def test_integrals():
    x = Symbol("x")
    for c in (Integral, Integral(x)):
        check(c)

#================== matrices ====================
from sympy.matrices.matrices import Matrix, SMatrix

def test_matrices():
    for c in (Matrix, Matrix([1,2,3]), SMatrix, SMatrix([[1,2],[3,4]])):
        check(c)

#================== ntheorie ====================
from sympy.ntheory.generate import Sieve

def test_ntheory():
    for c in (Sieve, Sieve()):
        check(c)

#================== physics =====================
from sympy.physics.paulialgebra import Pauli
from sympy.physics.units import Unit

def test_physics():
    for c in (Unit, Unit("meter", "m"), Pauli, Pauli(1)):
        check(c)

#================== plotting ====================
# XXX: These tests are not complete.

# these depend on ctypes, that are not in python2.4 by default, so XFAIled

@XFAIL
def test_plotting():
    from sympy.plotting.color_scheme import ColorGradient, ColorScheme
    from sympy.plotting.managed_window import ManagedWindow
    from sympy.plotting.plot import Plot, ScreenShot
    from sympy.plotting.plot_axes import PlotAxes, PlotAxesBase, PlotAxesFrame, PlotAxesOrdinate
    from sympy.plotting.plot_camera import PlotCamera
    from sympy.plotting.plot_controller import PlotController
    from sympy.plotting.plot_curve import PlotCurve
    from sympy.plotting.plot_interval import PlotInterval
    from sympy.plotting.plot_mode import PlotMode
    from sympy.plotting.plot_modes import Cartesian2D, Cartesian3D, Cylindrical,\
        ParametricCurve2D, ParametricCurve3D, ParametricSurface, Polar, Spherical
    from sympy.plotting.plot_object import PlotObject
    from sympy.plotting.plot_surface import PlotSurface
    from sympy.plotting.plot_window import PlotWindow
    for c in (ColorGradient, ColorGradient(0.2,0.4), ColorScheme, ManagedWindow,
              ManagedWindow, Plot, ScreenShot, PlotAxes, PlotAxesBase,
              PlotAxesFrame, PlotAxesOrdinate, PlotCamera, PlotController,
              PlotCurve, PlotInterval, PlotMode, Cartesian2D, Cartesian3D,
              Cylindrical, ParametricCurve2D, ParametricCurve3D,
              ParametricSurface, Polar, Spherical, PlotObject, PlotSurface,
              PlotWindow):
        check(c)

@XFAIL
def test_plotting2():
    from sympy.plotting.color_scheme import ColorGradient, ColorScheme
    from sympy.plotting.managed_window import ManagedWindow
    from sympy.plotting.plot import Plot, ScreenShot
    from sympy.plotting.plot_axes import PlotAxes, PlotAxesBase, PlotAxesFrame, PlotAxesOrdinate
    from sympy.plotting.plot_camera import PlotCamera
    from sympy.plotting.plot_controller import PlotController
    from sympy.plotting.plot_curve import PlotCurve
    from sympy.plotting.plot_interval import PlotInterval
    from sympy.plotting.plot_mode import PlotMode
    from sympy.plotting.plot_modes import Cartesian2D, Cartesian3D, Cylindrical,\
        ParametricCurve2D, ParametricCurve3D, ParametricSurface, Polar, Spherical
    from sympy.plotting.plot_object import PlotObject
    from sympy.plotting.plot_surface import PlotSurface
    from sympy.plotting.plot_window import PlotWindow
    check(ColorScheme("rainbow"))
    check(Plot(1,visible=False))
    check(PlotAxes())

#================== polys =======================
from sympy.polys.polynomial import IntegerPoly, Poly
from sympy.polys.rootfinding import RootOf, RootsOf, RootSum

def test_polys():
    x = Symbol("x")
    f = Poly(x, x)
    g = lambda x: x

    for c in (IntegerPoly, IntegerPoly(x, x), Poly, Poly(x, x)):
        check(c)

    for c in (RootOf, RootOf(f, 0), RootsOf, RootsOf(x, x), RootSum, RootSum(g, f)):
        check(c)

#================== printing ====================
from sympy.printing.latex import LatexPrinter
from sympy.printing.mathml import MathMLPrinter
from sympy.printing.pretty.pretty import PrettyPrinter
from sympy.printing.pretty.stringpict import prettyForm, stringPict
from sympy.printing.printer import Printer
from sympy.printing.python import PythonPrinter

def test_printing():
    for c in (LatexPrinter, LatexPrinter(), MathMLPrinter,
              PrettyPrinter, prettyForm, stringPict, stringPict("a"),
              Printer, Printer(), PythonPrinter, PythonPrinter()):
        check(c)

@XFAIL
def test_printing1():
    check(MathMLPrinter())

@XFAIL
def test_printing2():
    check(PrettyPrinter())

#================== series ======================
from sympy.series.gruntz import Limit2
from sympy.series.limits import Limit
from sympy.series.order import Order

def test_series():
    e = Symbol("e")
    x = Symbol("x")
    for c in (Limit2, Limit2(e, x, 1), Limit, Limit(e, x, 1), Order, Order(e)):
        check(c)

#================== statistics ==================
from sympy.statistics.distributions import ContinuousProbability, Normal, Sample, Uniform

def test_statistics():
    x = Symbol("x")
    y = Symbol("y")
    for c in (ContinuousProbability, ContinuousProbability(), Normal,
              Normal(x,y), Sample, Sample([1,3,4]), Uniform, Uniform(x,y)):
        check(c)

#================== concrete ==================
from sympy.concrete.products import Product
from sympy.concrete.summations import Sum
from sympy.concrete.sums_products import Sum2, _BigOperator

def test_concrete():
    x = Symbol("x")
    for c in (Product, Product(1,2), Sum, Sum(1), Sum2, Sum2(x,(x,2,4)),
              _BigOperator):
        check(c)
