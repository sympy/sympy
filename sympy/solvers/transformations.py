from __future__ import print_function, division

from sympy.core.function import (Lambda, expand, expand_complex)
from sympy.core.relational import Eq
                             arg, Piecewise, piecewise_fold)
from sympy.functions.elementary.trigonometric import (TrigonometricFunction,
                                                      HyperbolicFunction)
from sympy.sets import (FiniteSet, EmptySet, imageset, Interval, Intersection,
                        Union, ConditionSet)


def TR_product(in_set):
    # XXX: assuming it is function of only one variable
    x = in_set.condition.variables[0]
    cond = in_set.condition.expr
    if isinstance(cond, Eq):
        h = cond.lhs - cond.rhs
        if h.is_Mul:
            f = h.args[0]
            g = h.args[1]
            if f.is_finite and g.is_finite:
                return Union(ConditionSet(Lambda(x, Eq(f, 0)), in_set.base_set),
                             ConditionSet(Lambda(x, Eq(g, 0)), in_set.base_set))
    return in_set
