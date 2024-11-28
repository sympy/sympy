#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sympy                     import Basic, sympify
from sympy.core.expr           import Expr
from sympy.core.numbers        import Integer
from sympy.core.symbol         import Symbol
from sympy.physics             import units as u

# `Source txt2sympy1 <https://docs.djangoproject.com/en/1.9/howto/custom-model-fields/>`_ and
# `Source txt2sympy2 <https://stackoverflow.com/questions/15895819/how-to-parse-and-simplify-a-string-like-3cm-%C2%B5s%C2%B2-4e-4-sqmiles-km-h2-treatin>`_
def str2sympy(sympy_str):
    """Convert string to Sympy object."""

    if isinstance(sympy_str, Basic) or sympy_str is None:
        return sympy_str
    subs = {}
    for k, v in u.__dict__.items():
        if (isinstance(v, Expr) and v.has(u.Unit)) or isinstance(v, Integer):
            subs[Symbol(k)] = v
    if sympify(sympy_str) is None:
        return None
    else:
        return sympify(sympy_str).subs(subs)
