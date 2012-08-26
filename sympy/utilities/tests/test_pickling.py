import copy
import pickle
import warnings
import sys
from sympy.utilities.pytest import XFAIL

from sympy.core.basic import Atom, Basic
from sympy.core.core import BasicMeta, BasicType, ClassRegistry
from sympy.core.singleton import SingletonRegistry
from sympy.core.symbol import Dummy, Symbol, Wild
from sympy.core.numbers import (E, I, pi, oo, zoo, nan, Integer, Number,
        NumberSymbol, Rational, Float)
from sympy.core.relational import (Equality, GreaterThan, LessThan, Relational,
        StrictGreaterThan, StrictLessThan, Unequality)
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.core.function import Derivative, Function, FunctionClass, Lambda,\
        WildFunction
from sympy.core.sets import Interval
from sympy.core.multidimensional import vectorize
from sympy.functions import exp
#from sympy.core.ast_parser import SymPyParser, SymPyTransformer

from sympy.core.compatibility import callable
from sympy.utilities.exceptions import SymPyDeprecationWarning

from sympy import symbols, S

excluded_attrs = set(['_assumptions', '_mhash'])

def check(a, check_attr=True):
    """ Check that pickling and copying round-trips.
    """
    # The below hasattr() check will warn about is_Real in Python 2.5, so
    # disable this to keep the tests clean
    warnings.filterwarnings("ignore", category=SymPyDeprecationWarning)
    protocols = [0, 1, 2, copy.copy, copy.deepcopy]
    # Python 2.x doesn't support the third pickling protocol
    if sys.version_info[0] > 2:
        protocols.extend([3])
    for protocol in protocols:
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

        if not check_attr:
            continue
        def c(a, b, d):
            for i in d:
                if not hasattr(a, i) or i in excluded_attrs:
                    continue
                attr = getattr(a, i)
                if not hasattr(attr, "__call__"):
                    assert hasattr(b,i), i
                    assert getattr(b,i) == attr
        c(a,b,d1)
        c(b,a,d2)

    warnings.filterwarnings("default", category=SymPyDeprecationWarning)

#================== core =========================

def test_core_basic():
    for c in (Atom, Atom(),
              Basic, Basic(),
              # XXX: dynamically created types are not picklable
              # BasicMeta, BasicMeta("test", (), {}),
              # BasicType, BasicType("test", (), {}),
              ClassRegistry, ClassRegistry(),
              SingletonRegistry, SingletonRegistry()):
        check(c)

def test_core_symbol():
    # make the Symbol a unique name that doesn't class with any other
    # testing variable in this file since after this test the symbol
    # having the same name will be cached as noncommutative
    for c in (Dummy, Dummy("x", commutative=False), Symbol,
            Symbol("_issue_3130", commutative=False), Wild, Wild("x")):
        check(c)

def test_core_numbers():
    for c in (Integer(2), Rational(2, 3), Float("1.2")):
        check(c)

def test_core_relational():
    x = Symbol("x")
    y = Symbol("y")
    for c in (Equality, Equality(x,y), GreaterThan, GreaterThan(x, y),
              LessThan, LessThan(x,y), Relational, Relational(x,y),
              StrictGreaterThan, StrictGreaterThan(x,y), StrictLessThan,
              StrictLessThan(x,y), Unequality, Unequality(x,y)):
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

def test_Singletons():
    protocols = [0, 1, 2]
    if sys.version_info[0] > 2:
        protocols.extend([3])
    copiers = [copy.copy, copy.deepcopy]
    copiers += [lambda x: pickle.loads(pickle.dumps(x, proto))
            for proto in protocols]

    for obj in (Integer(-1), Integer(0), Integer(1), Rational(1, 2), pi, E, I,
            oo, -oo, zoo, nan, S.GoldenRatio, S.EulerGamma, S.Catalan,
            S.EmptySet, S.IdentityFunction):
        for func in copiers:
            assert func(obj) is obj


#================== functions ===================
from sympy.functions import (Piecewise, lowergamma, acosh,
        chebyshevu, chebyshevt, ln, chebyshevt_root, binomial, legendre,
        Heaviside, factorial, bernoulli, coth, tanh, assoc_legendre, sign,
        arg, asin, DiracDelta, re, rf, Abs, uppergamma, binomial, sinh, Ylm,
        cos, cot, acos, acot, gamma, bell, hermite, harmonic,
        LambertW, zeta, log, factorial, asinh, acoth, Zlm,
        cosh, dirichlet_eta, Eijk, loggamma, erf, ceiling, im, fibonacci,
        conjugate, tan, chebyshevu_root, floor, atanh, sqrt,
        RisingFactorial, sin, atan, ff, FallingFactorial, lucas, atan2,
        polygamma, exp)

def test_functions():
    one_var = (acosh, ln, Heaviside, factorial, bernoulli, coth, tanh,
            sign, arg, asin, DiracDelta, re, Abs, sinh, cos, cot, acos, acot,
            gamma, bell, harmonic, LambertW, zeta, log, factorial, asinh,
            acoth, cosh, dirichlet_eta, loggamma, erf, ceiling, im, fibonacci,
            conjugate, tan, floor, atanh, sin, atan, lucas, exp)
    two_var = (rf, ff, lowergamma, chebyshevu, chebyshevt, binomial,
            atan2, polygamma, hermite, legendre, uppergamma)
    x, y, z = symbols("x,y,z")
    others = (chebyshevt_root, chebyshevu_root, Eijk(x, y, z),
            Piecewise( (0, x<-1), (x**2, x<=1), (x**3, True)),
            assoc_legendre)
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
              Triangle, Triangle(p1,p2,p3)):
        check(c, check_attr = False)

#================== integrals ====================
from sympy.integrals.integrals import Integral

def test_integrals():
    x = Symbol("x")
    for c in (Integral, Integral(x)):
        check(c)

#==================== logic =====================
from sympy.core.logic import Logic

def test_logic():
    for c in (Logic, Logic(1)):
        check(c)

#================== matrices ====================
from sympy.matrices.matrices import Matrix, SparseMatrix

def test_matrices():
    for c in (Matrix, Matrix([1,2,3]), SparseMatrix, SparseMatrix([[1,2],[3,4]])):
        check(c)

#================== ntheory =====================
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
# XXX: These tests are not complete, so XFAIL them

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
from sympy.polys.polytools import Poly
from sympy.polys.polyclasses import DMP, DMF, ANP
from sympy.polys.rootoftools import RootOf, RootSum

from sympy.polys.domains import (
    PythonIntegerRing,
    SymPyIntegerRing,
    SymPyRationalField,
    PolynomialRing,
    FractionField,
    ExpressionDomain,
)

def test_polys():
    x = Symbol("X")

    ZZ = PythonIntegerRing()
    QQ = SymPyRationalField()

    for c in (Poly, Poly(x, x)):
        check(c)

    for c in (DMP, DMP([[ZZ(1)],[ZZ(2)],[ZZ(3)]], ZZ)):
        check(c)
    for c in (DMF, DMF(([ZZ(1),ZZ(2)], [ZZ(1),ZZ(3)]), ZZ)):
        check(c)
    for c in (ANP, ANP([QQ(1),QQ(2)], [QQ(1),QQ(2),QQ(3)], QQ)):
        check(c)

    for c in (PythonIntegerRing, PythonIntegerRing()):
        check(c)
    for c in (SymPyIntegerRing, SymPyIntegerRing()):
        check(c)
    for c in (SymPyRationalField, SymPyRationalField()):
        check(c)

    for c in (PolynomialRing, PolynomialRing(ZZ, 'x', 'y')):
        check(c)
    for c in (FractionField, FractionField(ZZ, 'x', 'y')):
        check(c)

    for c in (ExpressionDomain, ExpressionDomain()):
        check(c)

    from sympy.polys.domains import PythonRationalField

    for c in (PythonRationalField, PythonRationalField()):
        check(c)

    from sympy.polys.domains import HAS_GMPY

    if HAS_GMPY:
        from sympy.polys.domains import GMPYIntegerRing, GMPYRationalField

        for c in (GMPYIntegerRing, GMPYIntegerRing()):
            check(c)
        for c in (GMPYRationalField, GMPYRationalField()):
            check(c)

    f = x**3 + x + 3
    g = exp

    for c in (RootOf, RootOf(f, 0), RootSum, RootSum(f, g)):
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
from sympy.series.limits import Limit
from sympy.series.order import Order

def test_series():
    e = Symbol("e")
    x = Symbol("x")
    for c in (Limit, Limit(e, x, 1), Order, Order(e)):
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

def test_concrete():
    x = Symbol("x")
    for c in (Product, Product(x, (x, 2, 4)), Sum, Sum(x, (x, 2, 4))):
        check(c)
