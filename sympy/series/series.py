from __future__ import print_function, division

from sympy.core.sympify import sympify


def series(expr, x=None, x0=0, n=6, dir="+"):
    """Series expansion of expr around point `x = x0`.

       Parameters
       ==========

       expr : Expression
              The expression whose series is to be expanded.

       x    : Symbol
              It is the variable of the expression to be calculated.

       x0   : Value
              The value around which ``x`` is calculated. Can be any value
              from ``-oo`` to ``oo``.

       n    : Value
              The number of terms upto which the series is to be expanded.

       dir  : String
              Optional (default: "+")
              The series-expansion can be bi-directional. If ``dir="+"``,
              then (x->x0+). If ``dir="-", then (x->x0-). For infinite
              ``x0`` (``oo`` or ``-oo``), the ``dir`` argument is determined
              from the direction of the infinity (i.e., ``dir="-"`` for
              ``oo``).

       Examples
       ========

       >>> from sympy import Symbol, series, tan, oo
       >>> from sympy.abc import x
       >>> f = tan(x)
       >>> series(f, x, 2, 6, "+")
       tan(2) + (1 + tan(2)**2)*(x - 2) + (x - 2)**2*(tan(2)**3 + tan(2)) +
       (x - 2)**3*(1/3 + 4*tan(2)**2/3 + tan(2)**4) + (x - 2)**4*(tan(2)**5 +
       5*tan(2)**3/3 + 2*tan(2)/3) + (x - 2)**5*(2/15 + 17*tan(2)**2/15 +
       2*tan(2)**4 + tan(2)**6) + O((x - 2)**6, (x, 2))

       >>> series(f, x, 2, 3, "-")
       tan(2) + (2 - x)*(-tan(2)**2 - 1) + (2 - x)**2*(tan(2)**3 + tan(2))
       + O((x - 2)**3, (x, 2))

       >>> series(f, x, 2, oo, "+")
       Traceback (most recent call last):
         File "/opt/python/3.7.1/lib/python3.7/doctest.py", line 1329, in __run
           compileflags, 1), test.globs)
         File "<doctest sympy.series.series.series[5]>", line 1, in <module>
           series(f, x, 2, oo, "+")
         File "/home/travis/virtualenv/python3.7.1/lib/python3.7/site-packages/sympy-1.5.dev0-py3.7.egg/sympy/core/expr.py", line 2721, in series
           s = self.subs(x, rep).series(x, x0=0, n=n, dir='+', logx=logx)
         File "/home/travis/virtualenv/python3.7.1/lib/python3.7/site-packages/sympy-1.5.dev0-py3.7.egg/sympy/core/expr.py", line 2731, in series
           rv = self.subs(x, xpos).series(xpos, x0, n, dir, logx=logx)
         File "/home/travis/virtualenv/python3.7.1/lib/python3.7/site-packages/sympy-1.5.dev0-py3.7.egg/sympy/core/expr.py", line 2738, in series
           s1 = self._eval_nseries(x, n=n, logx=logx)
         File "/home/travis/virtualenv/python3.7.1/lib/python3.7/site-packages/sympy-1.5.dev0-py3.7.egg/sympy/functions/elementary/trigonometric.py", line 1094, in _eval_nseries
           return Function._eval_nseries(self, x, n=n, logx=logx)
         File "/home/travis/virtualenv/python3.7.1/lib/python3.7/site-packages/sympy-1.5.dev0-py3.7.egg/sympy/core/function.py", line 707, in _eval_nseries
           for i in range(n - 1):
       TypeError: 'Infinity' object cannot be interpreted as an integer

       Returns
       =======

       Series expansion of the expression about x0

       See Also
       ========

       See the docstring of Expr.series() for complete details of this wrapper.
    """
    expr = sympify(expr)
    return expr.series(x, x0, n, dir)
