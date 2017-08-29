# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from itertools import product
from sympy import Add, Symbol, sin, Abs, oo, Interval
from sympy.calculus.singularities import is_increasing, is_decreasing
from sympy.codegen.rewriting import Optimization

"""
This module collects classes useful for approimate rewriting of expressions.
This can be beneficial when generating numeric code for which performance is
of greater importance than precision (e.g. for preconditioners used in iterative
methods).
"""

class SumApprox(Optimization):
    """ Approximates sum by neglecting small terms

    If terms are expressions which can be determined to be monotonic, then
    bounds for those expressions are added.

    Parameters
    ----------
    bounds : dict
        Mapping expressions to length 2 tuple of bounds (low, high).
    reltol : number
        Threshold for when to ignore a term. Taken relative to the largest
        lower bound among bounds.

    Examples
    --------
    >>> from sympy import exp
    >>> from sympy.abc import x, y, z
    >>> from sympy.codegen.rewriting import optimize
    >>> from sympy.codegen.approximations import SumApprox
    >>> bounds = {x: (-1, 1), y: (1000, 2000), z: (-10, 3)}
    >>> sum_approx3 = SumApprox(bounds, reltol=1e-3)
    >>> sum_approx2 = SumApprox(bounds, reltol=1e-2)
    >>> sum_approx1 = SumApprox(bounds, reltol=1e-1)
    >>> expr = 3*(x + y + exp(z))
    >>> optimize(expr, [sum_approx3])
    3*(x + y + exp(z))
    >>> optimize(expr, [sum_approx2])
    3*y + 3*exp(z)
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
        for term in add.args:
            if term.is_number or term in self.bounds or len(term.free_symbols) != 1:
                continue
            fs, = term.free_symbols
            if fs not in self.bounds:
                continue
            intrvl = Interval(*self.bounds[fs])
            if is_increasing(term, intrvl, fs):
                self.bounds[term] = (
                    term.subs({fs: self.bounds[fs][0]}),
                    term.subs({fs: self.bounds[fs][1]})
                )
            elif is_decreasing(term, intrvl, fs):
                self.bounds[term] = (
                    term.subs({fs: self.bounds[fs][1]}),
                    term.subs({fs: self.bounds[fs][0]})
                )
            else:
                return add

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
