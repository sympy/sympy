# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from itertools import product
from sympy import Add, Symbol, sin, Abs, oo
from sympy.codegen.rewriting import Optimization

"""
This module collects classes useful for approimate rewriting of expressions.
This can be beneficial when generating numeric code for which performance is
of greater importance than precision (e.g. for preconditioners used in iterative
methods).
"""

class SumApprox(Optimization):
    """ Approximates sum by neglecting small terms

    Parameters
    ----------
    bounds : dict
        Mapping expressions to length 2 tuple of bounds (low, high).
    reltol : number
        Threshold for when to ignore a term. Taken relative to the largest
        lower bound among bounds.

    Examples
    --------
    >>> from sympy.abc import x, y, z
    >>> from sympy.codegen.rewriting import optimize
    >>> from sympy.codegen.approximations import SumApprox
    >>> bounds = {x: (-1, 1), y: (1000, 2000), z: (0, 20)}
    >>> sum_approx3 = SumApprox(bounds, reltol=1e-3)
    >>> sum_approx2 = SumApprox(bounds, reltol=1e-2)
    >>> sum_approx1 = SumApprox(bounds, reltol=1e-1)
    >>> expr = 3*(x + y + z)
    >>> optimize(expr, [sum_approx3])
    3*(x + y + z)
    >>> optimize(expr, [sum_approx2])
    3*y + 3*z
    >>> optimize(expr, [sum_approx1])
    3*y

    """

    def __init__(self, bounds, reltol, **kwargs):
        super(SumApprox, self).__init__(**kwargs)
        self.bounds = bounds
        self.reltol = reltol

    def __call__(self, expr):
        return expr.factor().replace(self.query, lambda arg: self.value(arg))

    def query(self, expr):
        return expr.is_Add

    def value(self, add):
        if all(term.is_number or term in self.bounds for term in add.args):
            bounds = [(term, term) if term.is_number else self.bounds[term] for term in add.args]
            largest_abs_guarantee = 0
            for lo, hi in bounds:
                if lo <= 0 <= hi:
                    continue
                largest_abs_guarantee = max(largest_abs_guarantee,
                                            min(abs(lo), abs(hi)))
            new_terms = []
            for term, (lo, hi) in zip(add.args, bounds):
                if max(abs(lo), abs(hi)) >= largest_abs_guarantee*self.reltol:
                    new_terms.append(term)
            return add.func(*new_terms)
        else:
            return add
