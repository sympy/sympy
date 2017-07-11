"""Core module. Provides the basic operations needed in sympy.
"""

from .add import Add
from .basic import Atom, Basic, preorder_traversal
from .cache import cacheit
from .containers import Dict, Tuple
from .evalf import N, PrecisionExhausted
from .evaluate import evaluate
from .expr import AtomicExpr, Expr, UnevaluatedExpr
from .exprtools import factor_nc, factor_terms, gcd_terms
from .function import Derivative, Function, FunctionClass, Lambda, PoleError, \
    Subs, WildFunction, count_ops, diff, expand, expand_complex, expand_func, \
    expand_log, expand_mul, expand_multinomial, expand_power_base, \
    expand_power_exp, expand_trig, nfloat
from .mod import Mod
from .mul import Mul, prod
from .multidimensional import vectorize
from .numbers import AlgebraicNumber, E, Float, I, Integer, Number, \
    NumberSymbol, Rational, RealNumber, comp, igcd, ilcm, mod_inverse, nan, \
    oo, pi, seterr, zoo
from .power import Pow, integer_nthroot
from .relational import Eq, Equality, Ge, GreaterThan, Gt, Le, LessThan, Lt, \
    Ne, Rel, StrictGreaterThan, StrictLessThan, Unequality
from .singleton import S
from .symbol import Dummy, Symbol, Wild, symbols, var
from .sympify import SympifyError, sympify

# expose singletons
Catalan = S.Catalan
EulerGamma = S.EulerGamma
GoldenRatio = S.GoldenRatio
