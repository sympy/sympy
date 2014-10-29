from __future__ import print_function, division

from sympy.core.add import Add
from sympy.core.singleton import S
from sympy.core.function import Function, ArgumentIndexError


def _eval_derivative_heurisch(f, s):
    # f(x).diff(s) -> x.diff(s) * f.fdiff(1)(s)
    i = 0
    l = []
    for a in f.args:
        i += 1
        da = a.diff(s)
        if da is S.Zero:
            continue
        try:
            df = f.heurisch_fdiff(i)
        except ArgumentIndexError:
            df = Function.fdiff(f, i)
        l.append(df * da)
    return Add(*l)

def replace_heurisch_diff(function):
    cls = function.func
    overrides = dict(_eval_derivative=_eval_derivative_heurisch,
                     _old_class=cls)
    newcls = type(cls.__name__, (cls,), overrides)
    return newcls(*function.args)

def restore_heurisch_diff(function):
    cls = function.func
    oldcls = cls._old_class
    return oldcls(*function.args)

def has_heurisch_diff(f):
    return hasattr(f, "heurisch_fdiff")
