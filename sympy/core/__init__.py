"""Core module. Provides the basic operations needed in sympy.
"""

__all__ = []

from .sympify import sympify, SympifyError
__all__ += ["sympify", "SympifyError"]

from .cache import cacheit
__all__ += ["cacheit"]

from .basic import Basic, Atom, preorder_traversal
__all__ += ["Basic", "Atom", "preorder_traversal"]

from .singleton import S
__all__ += ["S"]

from .expr import Expr, AtomicExpr, UnevaluatedExpr
__all__ += ["Expr", "AtomicExpr", "UnevaluatedExpr"]

from .symbol import Symbol, Wild, Dummy, symbols, var
__all__ += ["Symbol", "Wild", "Dummy", "symbols", "var"]

from .numbers import (
    Number, Float, Rational, Integer, NumberSymbol,
    RealNumber, igcd, ilcm, seterr, E, I, nan, oo, pi, zoo,
    AlgebraicNumber, comp, mod_inverse
)
__all__ += [
    "Number", "Float", "Rational", "Integer", "NumberSymbol",
    "RealNumber", "igcd", "ilcm", "seterr", "E", "I", "nan", "oo", "pi", "zoo",
    "AlgebraicNumber", "comp", "mod_inverse"
]

from .power import Pow, integer_nthroot
__all__ += ["Pow", "integer_nthroot"]

from .mul import Mul, prod
__all__ += ["Mul", "prod"]

from .add import Add
__all__ += ["Add"]

from .mod import Mod
__all__ += ["Mod"]

from .relational import (
    Rel, Eq, Ne, Lt, Le, Gt, Ge,
    Equality, GreaterThan, LessThan, Unequality,
    StrictGreaterThan, StrictLessThan
)
__all__ += [
    "Rel", "Eq", "Ne", "Lt", "Le", "Gt", "Ge",
    "Equality", "GreaterThan", "LessThan", "Unequality",
    "StrictGreaterThan", "StrictLessThan"
]

from .multidimensional import vectorize
__all__ += ["vectorize"]

from .function import (
    Lambda, WildFunction, Derivative, diff, FunctionClass,
    Function, Subs, expand, PoleError, count_ops,
    expand_mul, expand_log, expand_func,
    expand_trig, expand_complex, expand_multinomial, nfloat,
    expand_power_base, expand_power_exp
)
__all__ += [
    "Lambda", "WildFunction", "Derivative", "diff", "FunctionClass",
    "Function", "Subs", "expand", "PoleError", "count_ops",
    "expand_mul", "expand_log", "expand_func",
    "expand_trig", "expand_complex", "expand_multinomial", "nfloat",
    "expand_power_base", "expand_power_exp"
]

from .evalf import PrecisionExhausted, N
__all__ += ["PrecisionExhausted", "N", "evalf"]

from .containers import Tuple, Dict
__all__ += ["Tuple", "Dict"]

from .exprtools import gcd_terms, factor_terms, factor_nc
__all__ += ["gcd_terms", "factor_terms", "factor_nc"]

from .evaluate import evaluate
__all__ += ["evaluate"]


# expose singletons
Catalan = S.Catalan
EulerGamma = S.EulerGamma
GoldenRatio = S.GoldenRatio
__all__ += ["Catalan", "EulerGamma", "GoldenRatio"]
