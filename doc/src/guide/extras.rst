
.. include:: ../definitions.def

======
Extras
======

Solutions
=========

    >>> from sympy import init_printing
    >>> init_printing(use_unicode=True, no_global=True)

Arithmetic operators
--------------------

.. _solution_arith_op_1:

Solution 1
~~~~~~~~~~

::

    >>> from sympy import Symbol, Add, Poly

    >>> x = Symbol('x')

    >>> Add(*[ x**i for i in range(0, 10+1) ])
     10    9    8    7    6    5    4    3    2
    x   + x  + x  + x  + x  + x  + x  + x  + x  + x + 1

    >>> sum([ x**i for i in range(0, 10+1) ])
     10    9    8    7    6    5    4    3    2
    x   + x  + x  + x  + x  + x  + x  + x  + x  + x + 1

    >>> Poly([1]*11, x).as_expr()
     10    9    8    7    6    5    4    3    2
    x   + x  + x  + x  + x  + x  + x  + x  + x  + x + 1

    >>> def build_poly_1(n):
    ...     return Add(*[ x**i for i in range(0, n+1) ])
    ...

    >>> build_poly_1(1)
    x + 1
    >>> build_poly_1(2)
     2
    x  + x + 1
    >>> build_poly_1(3)
     3    2
    x  + x  + x + 1

    >>> def build_poly_2(n):
    ...     if n > 0:
    ...         return Add(*[ x**i for i in range(0, n+1) ])
    ...     else:
    ...         return Add(*[ x**i for i in range(0, n-1, -1) ])
    ...

    >>> build_poly_2(1)
    x + 1
    >>> build_poly_2(0)
    1
    >>> build_poly_2(-1)
        1
    1 + ─
        x

.. _solution_arith_op_2:

Solution 2
~~~~~~~~~~

::

    >>> from sympy import Symbol

    >>> x = Symbol('x')

    >>> def nested_power(expr, n):
    ...     if not n:
    ...         return expr
    ...     else:
    ...         return expr**nested_power(expr, n-1)
    ...

    >>> nested_power(x, 1)
     x
    x
    >>> nested_power(x, 2)
     ⎛ x⎞
     ⎝x ⎠
    x
    >>> nested_power(x, 3)
     ⎛ ⎛ x⎞⎞
     ⎜ ⎝x ⎠⎟
     ⎝x    ⎠
    x

Building blocks of expressions
------------------------------

.. _solution_blocks_1:

Solution 1
~~~~~~~~~~

::

    >>> from sympy import var, exp, sin, Integral, Eq, pprint

    >>> var('a,b,c,x,n,C')
    (a, b, c, x, n, C)

    >>> expressions = [x**n, a*x**2 + b*x + c, 1/x, 1/(x**2 + 1), exp(a*x), sin(a*x) + b]
    >>> expressions
    ⎡ n     2            1    1      a⋅x              ⎤
    ⎢x , a⋅x  + b⋅x + c, ─, ──────, ℯ   , b + sin(a⋅x)⎥
    ⎢                    x   2                        ⎥
    ⎣                       x  + 1                    ⎦

    >>> integrals_table = []

    >>> for expr in expressions:
    ...     integral = Integral(expr, x)
    ...     integrals_table.append((integral, integral.doit() + C))
    ...

    >>> for integral, antiderivative in integrals_table:
    ...     pprint(Eq(integral, antiderivative), use_unicode=True)
    ...
    ⌠              n + 1
    ⎮  n          x
    ⎮ x  dx = C + ──────
    ⌡             n + 1
    ⌠                              3      2
    ⎮ ⎛   2          ⎞          a⋅x    b⋅x
    ⎮ ⎝a⋅x  + b⋅x + c⎠ dx = C + ──── + ──── + c⋅x
    ⌡                            3      2
    ⌠
    ⎮ 1
    ⎮ ─ dx = C + log(x)
    ⎮ x
    ⌡
    ⌠
    ⎮   1
    ⎮ ────── dx = C + atan(x)
    ⎮  2
    ⎮ x  + 1
    ⌡
    ⌠                a⋅x
    ⎮  a⋅x          ℯ
    ⎮ ℯ    dx = C + ────
    ⌡                a
    ⌠                               cos(a⋅x)
    ⎮ (b + sin(a⋅x)) dx = C + b⋅x - ────────
    ⌡                                  a

Foreign types in SymPy
----------------------

.. _solution_foreign_1:

Solution 1
~~~~~~~~~~

::

    >>> from sympy.core.sympify import sympify, converter # doctest: +SKIP
    >>> from sympy import Rational # doctest: +SKIP
    >>> from gmpy import mpq # doctest: +SKIP

    >>> def mpq_to_Rational(obj):
    ...     return Rational(obj.numer(), obj.denom())
    ...

    >>> converter[type(mpq(1))] = mpq_to_Rational # doctest: +SKIP

    >>> sympify(mpq(1, 2)) # doctest: +SKIP
    1/2
    >>> type(_) # doctest: +SKIP
    <class 'sympy.core.numbers.Half'>

.. _solution_foreign_2:

Solution 2
~~~~~~~~~~

::

    >>> from sympy.core.sympify import converter, sympify, SympifyError # doctest: +SKIP
    >>> from sympy import Tuple # doctest: +SKIP
    >>> from numpy import array, ndarray # doctest: +SKIP

    >>> def ndarray_to_Tuple(obj):
    ...     if len(obj.shape) == 1:
    ...         return Tuple(*obj)
    ...     else:
    ...         raise SympifyError("only row NumPy arrays are allowed")
    ...

    >>> converter[ndarray] = ndarray_to_Tuple # doctest: +SKIP

    >>> sympify(array([1, 2, 3])) # doctest: +SKIP
    (1, 2, 3)

    >>> sympify(array([[1], [2], [3]])) # doctest: +SKIP
    Traceback (most recent call last):
    ...
    SympifyError: SympifyError: 'only row NumPy arrays are allowed'

The role of symbols
-------------------

.. _solution_symbols_1:

Solution 1
~~~~~~~~~~

::

    >>> from sympy import Add, Symbol, symbols, numbered_symbols

    >>> def build_expression_1(name, n):
    ...     return Add(*[ Symbol('%s_%d' % (name, i))**i for i in range(1, n+1) ])
    ...

    >>> def build_expression_2(name, n):
    ...     X = symbols('%s1:%d' % (name, n+1))
    ...     return Add(*[ x**i for x, i in zip(X, range(1, n+1)) ])
    ...

    >>> def build_expression_3(name, n):
    ...     X = numbered_symbols(name, start=1)
    ...     return Add(*[ x**i for x, i in zip(X, range(1, n+1)) ])
    ...

    >>> build_expression_1('x', 5)
           2     3     4     5
    x₁ + x₂  + x₃  + x₄  + x₅
    >>> build_expression_2('y', 5)
           2     3     4     5
    y₁ + y₂  + y₃  + y₄  + y₅
    >>> build_expression_2('z', 5)
           2     3     4     5
    z₁ + z₂  + z₃  + z₄  + z₅

Obtaining parts of expressions
------------------------------

.. _solution_parts_1:

Solution 1
~~~~~~~~~~

::

    >>> from sympy.core import Atom, sympify

    >>> def depth_1(expr):
    ...     expr = sympify(expr)
    ...
    ...     if isinstance(expr, Atom):
    ...         return 1
    ...     else:
    ...         return 1 + max([ depth_1(arg) for arg in expr.args ])
    ...

    >>> depth_1(117)
    1

    >>> def depth_2(expr):
    ...     def _depth(expr):
    ...         if isinstance(expr, Atom):
    ...             return 1
    ...         else:
    ...             return 1 + max([ _depth(arg) for arg in expr.args ])
    ...
    ...     return _depth(sympify(expr))
    ...

    >>> depth_2(117)
    1

.. _solution_parts_2:

Solution 2
~~~~~~~~~~

::

    >>> from sympy import Atom, Basic, sin, cos
    >>> var('x,y')
    (x, y)

    >>> def depth(expr):
    ...     if isinstance(expr, Atom):
    ...         return 1
    ...     else:
    ...         if isinstance(expr, Basic):
    ...             args = expr.args
    ...         elif hasattr(expr, '__iter__'):
    ...             args = expr
    ...         else:
    ...             raise ValueError("can't compute the depth of %s" % expr)
    ...
    ...         return 1 + max([ depth(arg) for arg in args ])
    ...

    >>> depth(x)
    1
    >>> depth(x + y)
    2
    >>> depth([x, y])
    2
    >>> depth((x, sin(x)))
    3
    >>> depth(set([(x, sin(x), sin(cos(x)))]))
    5

Immutability of expressions
---------------------------

.. _solution_immutability_1:

Solution 1
~~~~~~~~~~

::

    >>> var('x')
    x

    >>> x**x
     x
    x
    >>> _.subs(x, x**x)
        ⎛ x⎞
        ⎝x ⎠
    ⎛ x⎞
    ⎝x ⎠
    >>> _.subs(x, x**x)
              ⎛    ⎛ x⎞⎞
              ⎜    ⎝x ⎠⎟
              ⎜⎛ x⎞    ⎟
              ⎝⎝x ⎠    ⎠
    ⎛    ⎛ x⎞⎞
    ⎜    ⎝x ⎠⎟
    ⎜⎛ x⎞    ⎟
    ⎝⎝x ⎠    ⎠
    >>> _.subs(x, x**x)
                          ⎛          ⎛    ⎛ x⎞⎞⎞
                          ⎜          ⎜    ⎝x ⎠⎟⎟
                          ⎜          ⎜⎛ x⎞    ⎟⎟
                          ⎜          ⎝⎝x ⎠    ⎠⎟
                          ⎜⎛    ⎛ x⎞⎞          ⎟
                          ⎜⎜    ⎝x ⎠⎟          ⎟
                          ⎜⎜⎛ x⎞    ⎟          ⎟
                          ⎝⎝⎝x ⎠    ⎠          ⎠
    ⎛          ⎛    ⎛ x⎞⎞⎞
    ⎜          ⎜    ⎝x ⎠⎟⎟
    ⎜          ⎜⎛ x⎞    ⎟⎟
    ⎜          ⎝⎝x ⎠    ⎠⎟
    ⎜⎛    ⎛ x⎞⎞          ⎟
    ⎜⎜    ⎝x ⎠⎟          ⎟
    ⎜⎜⎛ x⎞    ⎟          ⎟
    ⎝⎝⎝x ⎠    ⎠          ⎠


    >>> expr = x**x

    >>> for i in range(3):
    ...     expr = (expr.base**expr.base)**(expr.exp**expr.exp)
    ...

    >>> expr
                          ⎛          ⎛    ⎛ x⎞⎞⎞
                          ⎜          ⎜    ⎝x ⎠⎟⎟
                          ⎜          ⎜⎛ x⎞    ⎟⎟
                          ⎜          ⎝⎝x ⎠    ⎠⎟
                          ⎜⎛    ⎛ x⎞⎞          ⎟
                          ⎜⎜    ⎝x ⎠⎟          ⎟
                          ⎜⎜⎛ x⎞    ⎟          ⎟
                          ⎝⎝⎝x ⎠    ⎠          ⎠
    ⎛          ⎛    ⎛ x⎞⎞⎞
    ⎜          ⎜    ⎝x ⎠⎟⎟
    ⎜          ⎜⎛ x⎞    ⎟⎟
    ⎜          ⎝⎝x ⎠    ⎠⎟
    ⎜⎛    ⎛ x⎞⎞          ⎟
    ⎜⎜    ⎝x ⎠⎟          ⎟
    ⎜⎜⎛ x⎞    ⎟          ⎟
    ⎝⎝⎝x ⎠    ⎠          ⎠


Turning strings into expressions
--------------------------------

.. _solution_sympify_1:

Solution 1
~~~~~~~~~~

.. ::

Customizing built-in printers
-----------------------------

.. _solution_custom_printers_1:

Solution 1
~~~~~~~~~~

.. TODO: adjust this solution

::

    >>> from sympy import Symbol, Lambda, sin

    >>> from sympy.printing.pretty.pretty import PrettyPrinter
    >>> from sympy.printing.pretty.stringpict import prettyForm

    >>> x = Symbol('x')

    >>> class LambdaPrettyPrinter(PrettyPrinter):
    ...     def _print_Lambda(self, expr):
    ...         args = expr.args[0].args + (expr.args[1],)
    ...         pform_head = prettyForm('Lambda')
    ...         pform_tail = self._print_seq(args, '(', ')')
    ...         pform = prettyForm(*pform_head.right(pform_tail))
    ...         return pform
    ...

    >>> LambdaPrettyPrinter().doprint(Lambda(x, sin(1/x)))
                 /1\
    Lambda(x, sin|-|)
                 \x/

.. _solution_custom_printers_2:

Solution 2
~~~~~~~~~~

::

    >>> from sympy.printing.pretty.pretty import PrettyPrinter
    >>> from sympy.printing.pretty.stringpict import prettyForm, stringPict

    >>> class PolyPrettyPrinter(PrettyPrinter):
    ...     def _print_Poly(self, poly):
    ...         expr = poly.as_expr()
    ...         gens = list(poly.gens)
    ...         domain = poly.get_domain()
    ...
    ...         pform = self._print(expr)
    ...         pform = prettyForm(*stringPict.next(pform, ", "))
    ...
    ...         for gen in gens:
    ...             pform = prettyForm(*stringPict.next(pform, self._print(gen)))
    ...             pform = prettyForm(*stringPict.next(pform, ", "))
    ...
    ...         pform = prettyForm(*stringPict.next(pform, "domain="))
    ...         pform = prettyForm(*stringPict.next(pform, self._print(domain)))
    ...
    ...         pform = prettyForm(*pform.parens("(", ")", ifascii_nougly=True))
    ...         pform = prettyForm(*prettyForm('Poly').right(pform))
    ...
    ...         return pform
    ...

    >>> PolyPrettyPrinter().doprint(Poly(x**2 + 1))
          2
    Poly(x  + 1, x, domain=ZZ)

    >>> PolyPrettyPrinter().doprint(Poly(x**2 + y/2))
          2   y
    Poly(x  + -, x, y, domain=QQ)
              2

Implementing printers from scratch
----------------------------------

.. _solution_new_printers_1:

Solution 1
~~~~~~~~~~

::

    >>> class ExtendedMathematicaPrinter(MathematicaPrinter): # doctest: +SKIP
    ...     def _print_Pi(self, expr):
    ...         return 'Pi'
    ...

    >>> MathematicaPrinter().doprint(pi) # doctest: +SKIP
    pi
    >>> ExtendedMathematicaPrinter().doprint(pi) # doctest: +SKIP
    Pi

.. _solution_new_printers_2:

Solution 2
~~~~~~~~~~

::

    >>> class ExtendedMathematicaPrinter(MathematicaPrinter): # doctest: +SKIP
    ...     def _print_Mul(self, expr):
    ...         prec = precedence(expr)
    ...         return "*".join([ self.parenthesize(arg, prec) for arg in expr.args ])
    ...

    >>> MathematicaPrinter().doprint(2*sin(x)) # doctest: +SKIP
    2*sin(x)
    >>> ExtendedMathematicaPrinter().doprint(2*sin(x)) # doctest: +SKIP
    2*Sin[x]

Code generation
---------------

.. _solution_codegen_1:

Solution 1
~~~~~~~~~~

::

    >>> from sympy.utilities.codegen import codegen
    >>> from sympy import Symbol, chebyshevt

    >>> x = Symbol('x')

    >>> print(codegen(("chebyshevt_20", chebyshevt(20, x)), "C", "file")[0][1])
    /******************************************************************************
     *                    Code generated with sympy 0.7.2-git                     *
     *                                                                            *
     *              See http://www.sympy.org/ for more information.               *
     *                                                                            *
     *                       This file is part of 'project'                       *
     ******************************************************************************/
    #include "file.h"
    #include <math.h>
    <BLANKLINE>
    double chebyshevt_20(double x) {
    <BLANKLINE>
       return 524288*pow(x, 20) - 2621440*pow(x, 18) + 5570560*pow(x, 16) - 6553600*pow(x, 14) + 4659200*pow(x, 12) - 2050048*pow(x, 10) + 549120*pow(x, 8) - 84480*pow(x, 6) + 6600*pow(x, 4) - 200*pow(x, 2) + 1;
    <BLANKLINE>
    }

.. _solution_codegen_2:

Solution 2
~~~~~~~~~~

::

    >>> from sympy.utilities.codegen import codegen
    >>> from sympy import Symbol, chebyshevt

    >>> x = Symbol('x')

    >>> functions = [ ("chebyshevt_%d" % i, chebyshevt(i, x)) for i in range(0, 10) ]
    >>> print(codegen(functions, "F95", "file")[0][1])
    !******************************************************************************
    !*                    Code generated with sympy 0.7.2-git                     *
    !*                                                                            *
    !*              See http://www.sympy.org/ for more information.               *
    !*                                                                            *
    !*                       This file is part of 'project'                       *
    !******************************************************************************
    <BLANKLINE>
    REAL*8 function chebyshevt_0()
    implicit none
    <BLANKLINE>
    chebyshevt_0 = 1
    <BLANKLINE>
    end function
    <BLANKLINE>
    REAL*8 function chebyshevt_1(x)
    implicit none
    REAL*8, intent(in) :: x
    <BLANKLINE>
    chebyshevt_1 = x
    <BLANKLINE>
    end function
    <BLANKLINE>
    REAL*8 function chebyshevt_2(x)
    implicit none
    REAL*8, intent(in) :: x
    <BLANKLINE>
    chebyshevt_2 = 2*x**2 - 1
    <BLANKLINE>
    end function
    <BLANKLINE>
    REAL*8 function chebyshevt_3(x)
    implicit none
    REAL*8, intent(in) :: x
    <BLANKLINE>
    chebyshevt_3 = 4*x**3 - 3*x
    <BLANKLINE>
    end function
    <BLANKLINE>
    REAL*8 function chebyshevt_4(x)
    implicit none
    REAL*8, intent(in) :: x
    <BLANKLINE>
    chebyshevt_4 = 8*x**4 - 8*x**2 + 1
    <BLANKLINE>
    end function
    <BLANKLINE>
    REAL*8 function chebyshevt_5(x)
    implicit none
    REAL*8, intent(in) :: x
    <BLANKLINE>
    chebyshevt_5 = 16*x**5 - 20*x**3 + 5*x
    <BLANKLINE>
    end function
    <BLANKLINE>
    REAL*8 function chebyshevt_6(x)
    implicit none
    REAL*8, intent(in) :: x
    <BLANKLINE>
    chebyshevt_6 = 32*x**6 - 48*x**4 + 18*x**2 - 1
    <BLANKLINE>
    end function
    <BLANKLINE>
    REAL*8 function chebyshevt_7(x)
    implicit none
    REAL*8, intent(in) :: x
    <BLANKLINE>
    chebyshevt_7 = 64*x**7 - 112*x**5 + 56*x**3 - 7*x
    <BLANKLINE>
    end function
    <BLANKLINE>
    REAL*8 function chebyshevt_8(x)
    implicit none
    REAL*8, intent(in) :: x
    <BLANKLINE>
    chebyshevt_8 = 128*x**8 - 256*x**6 + 160*x**4 - 32*x**2 + 1
    <BLANKLINE>
    end function
    <BLANKLINE>
    REAL*8 function chebyshevt_9(x)
    implicit none
    REAL*8, intent(in) :: x
    <BLANKLINE>
    chebyshevt_9 = 256*x**9 - 576*x**7 + 432*x**5 - 120*x**3 + 9*x
    <BLANKLINE>
    end function

Partial fraction decomposition
------------------------------

.. _solution_partfrac_1:

Solution 1
~~~~~~~~~~

* `\frac{3 x + 5}{(2 x + 1)^2}`::

    >>> from sympy import factor, together, apart, numer, denom, cancel
    >>> from sympy.solvers import solve_undetermined_coeffs

    >>> f = (3*x + 5)/(2*x + 1)**2; f
     3⋅x + 5
    ──────────
             2
    (2⋅x + 1)
    >>> var('A:B')
    (A, B)
    >>> p1 = A/(2*x + 1)
    >>> p2 = B/(2*x + 1)**2
    >>> p1, p2
    ⎛   A         B     ⎞
    ⎜───────, ──────────⎟
    ⎜2⋅x + 1           2⎟
    ⎝         (2⋅x + 1) ⎠
    >>> h = sum((p1, p2));h
       A          B
    ─────── + ──────────
    2⋅x + 1            2
              (2⋅x + 1)
    >>> hcombined = factor(together(h)); hcombined
    2⋅A⋅x + A + B
    ─────────────
               2
      (2⋅x + 1)
    >>> eq = Eq(hcombined, f);eq
    2⋅A⋅x + A + B    3⋅x + 5
    ───────────── = ──────────
               2             2
      (2⋅x + 1)     (2⋅x + 1)
    >>> coeffeq = Eq(numer(eq.lhs), numer(eq.rhs)); coeffeq
    2⋅A⋅x + A + B = 3⋅x + 5
    >>> sol = solve_undetermined_coeffs(coeffeq, [A, B], x); sol
    {A: 3/2, B: 7/2}
    >>> solution = h.subs(sol); solution
         3             7
    ─────────── + ────────────
    2⋅(2⋅x + 1)              2
                  2⋅(2⋅x + 1)
    >>> apart(f, x)
         3             7
    ─────────── + ────────────
    2⋅(2⋅x + 1)              2
                  2⋅(2⋅x + 1)
    >>> solution == apart(f, x)
    True

* `\frac{3 x + 5}{(u x + v)^2}`::

    >>> var('x,u,v')
    (x, u, v)
    >>> f = (3*x + 5)/(u*x + v)**2; f
     3⋅x + 5
    ──────────
             2
    (u⋅x + v)
    >>> var('A:B')
    (A, B)
    >>> p1 = A/(u*x + v)
    >>> p2 = B/(u*x + v)**2
    >>> p1, p2
    ⎛   A         B     ⎞
    ⎜───────, ──────────⎟
    ⎜u⋅x + v           2⎟
    ⎝         (u⋅x + v) ⎠
    >>> h = sum((p1, p2));h
       A          B
    ─────── + ──────────
    u⋅x + v            2
              (u⋅x + v)
    >>> hcombined = factor(together(h)); hcombined
    A⋅u⋅x + A⋅v + B
    ───────────────
                2
       (u⋅x + v)
    >>> eq = Eq(hcombined, f); eq
    A⋅u⋅x + A⋅v + B    3⋅x + 5
    ─────────────── = ──────────
                2              2
       (u⋅x + v)      (u⋅x + v)
    >>> coeffeq = Eq(numer(eq.lhs), numer(eq.rhs)); coeffeq
    A⋅u⋅x + A⋅v + B = 3⋅x + 5
    >>> lhs, rhs = Poly(coeffeq.lhs, x), Poly(coeffeq.rhs, x)
    >>> lhs
    Poly(A*u*x + A*v + B, x, domain='ZZ[u,v,A,B]')
    >>> rhs
    Poly(3*x + 5, x, domain='ZZ')
    >>> coeffeqs = [ Eq(lhs.nth(i), rhs.nth(i)) for i in range(2) ]; coeffeqs
    [A⋅v + B = 5, A⋅u = 3]
    >>> sol = solve_undetermined_coeffs(coeffeq, [A, B], x)
    >>> sol[B] = cancel(sol[B])
    >>> sol
    ⎧   3     5⋅u - 3⋅v⎫
    ⎨A: ─, B: ─────────⎬
    ⎩   u         u    ⎭
    >>> solution = h.subs(sol); solution
     5⋅u - 3⋅v          3
    ──────────── + ───────────
               2   u⋅(u⋅x + v)
    u⋅(u⋅x + v)
    >>> apart(f, x) # Note, you must include x!
     5⋅u - 3⋅v          3
    ──────────── + ───────────
               2   u⋅(u⋅x + v)
    u⋅(u⋅x + v)
    >>> solution == apart(f, x)
    True

* `\frac{(3 x + 5)^2}{(2 x + 1)^2}`::

    >>> from sympy import div, collect, solve

    >>> f = (3*x + 5)**2/(2*x + 1)**2;f
             2
    (3⋅x + 5)
    ──────────
             2
    (2⋅x + 1)
    >>> # Here, deg(p) == deg(q), so we have to use div()
    >>> p, q = numer(f), denom(f)
    >>> d, r = div(p, q) # p = d*q + r, i.e., f = p/q = d + r/q
    >>> d, r
    (9/4, 21⋅x + 91/4)
    >>> solution = d # This is the polynomial part of the solution
    >>> var('A:B')
    (A, B)
    >>> p1 = A/(2*x + 1)
    >>> p2 = B/(2*x + 1)**2
    >>> p1, p2
    ⎛   A         B     ⎞
    ⎜───────, ──────────⎟
    ⎜2⋅x + 1           2⎟
    ⎝         (2⋅x + 1) ⎠
    >>> h = sum((p1, p2));h
       A          B
    ─────── + ──────────
    2⋅x + 1            2
              (2⋅x + 1)
    >>> hcombined = factor(together(h));hcombined
    2⋅A⋅x + A + B
    ─────────────
               2
      (2⋅x + 1)
    >>> eq = Eq(hcombined, r/q);eq
    2⋅A⋅x + A + B   21⋅x + 91/4
    ───────────── = ───────────
               2              2
      (2⋅x + 1)      (2⋅x + 1)
    >>> coeffeq = Eq(numer(eq.lhs), r); coeffeq # No need to call numer() here; we know it's r
    2⋅A⋅x + A + B = 21⋅x + 91/4
    >>> coeffeqs = collect(coeffeq.lhs - coeffeq.rhs, x, evaluate=False)
    >>> # Note, since we're just subtracting coeffeq.lhs from coeffeq.rhs, we could have
    >>> # done it without using Eq() in the first place.
    >>> sol = solve(coeffeqs.values())
    >>> solution += h.subs(sol);solution # Note the +=
    9        21            49
    ─ + ─────────── + ────────────
    4   2⋅(2⋅x + 1)              2
                      4⋅(2⋅x + 1)
    >>> apart(f, x)
    9        21            49
    ─ + ─────────── + ────────────
    4   2⋅(2⋅x + 1)              2
                      4⋅(2⋅x + 1)
    >>> solution == apart(f, x)
    True

.. _solution_partfrac_2:

Solution 2
~~~~~~~~~~

Even though you can use :func:`Expr.coeff` to get the coefficient of
`x^n` in a polynomial for `n > 0`, it will not work for `n = 0`. This
is because :func:`Expr.coeff` considers any term with no numerical
coefficient to be the coefficient of 1. From the docstring::

    You can select terms with no rational coefficient:
    >>> (x+2*y).coeff(1)
    x
    >>> (3+2*x+4*x**2).coeff(1)
    0

Deriving trigonometric identities
---------------------------------

.. _solution_trig_1:

Solution 1
~~~~~~~~~~

::

    >>> var('a,b')
    (a, b)
    >>> f = sin(a + b);f
    sin(a + b)
    >>> g = f.series(a, 0, 10)
    >>> g
                         2           3           4           5           6
                        a ⋅sin(b)   a ⋅cos(b)   a ⋅sin(b)   a ⋅cos(b)   a ⋅sin(b)
    sin(b) + a⋅cos(b) - ───────── - ───────── + ───────── + ───────── - ─────────
                            2           6           24         120         720
    <BLANKLINE>
       7           8           9
      a ⋅cos(b)   a ⋅sin(b)   a ⋅cos(b)    ⎛ 10⎞
    - ───────── + ───────── + ───────── + O⎝a  ⎠
         5040       40320       362880
    >>> g = collect(g.removeO(), [sin(b), cos(b)])
    >>> g
    ⎛   8      6    4    2    ⎞          ⎛   9       7      5    3    ⎞
    ⎜  a      a    a    a     ⎟          ⎜  a       a      a    a     ⎟
    ⎜───── - ─── + ── - ── + 1⎟⋅sin(b) + ⎜────── - ──── + ─── - ── + a⎟⋅cos(b)
    ⎝40320   720   24   2     ⎠          ⎝362880   5040   120   6     ⎠
    >>> sina = sin(a).series(a, 0, 10).removeO()
    >>> sina
       9       7      5    3
      a       a      a    a
    ────── - ──── + ─── - ── + a
    362880   5040   120   6
    >>> h = g.subs(sina, sin(a))
    >>> h
    ⎛   8      6    4    2    ⎞
    ⎜  a      a    a    a     ⎟
    ⎜───── - ─── + ── - ── + 1⎟⋅sin(b) + sin(a)⋅cos(b)
    ⎝40320   720   24   2     ⎠
    >>> cosa = cos(a).series(a, 0, 10).removeO()
    >>> cosa
       8      6    4    2
      a      a    a    a
    ───── - ─── + ── - ── + 1
    40320   720   24   2
    >>> h = h.subs(cosa, cos(a))
    >>> h
    sin(a)⋅cos(b) + sin(b)⋅cos(a)
    >>> eq = Eq(f, h)
    >>> eq
    sin(a + b) = sin(a)⋅cos(b) + sin(b)⋅cos(a)
    >>> eq.lhs.expand(trig=True) == eq.rhs
    True

.. _solution_trig_2:

Solution 2
~~~~~~~~~~

::

    >>> var('a,b')
    (a, b)
    >>> f = cos(a + b);f
    cos(a + b)
    >>> g = f.series(a, 0, 10)
    >>> g
                         2           3           4           5           6
                        a ⋅cos(b)   a ⋅sin(b)   a ⋅cos(b)   a ⋅sin(b)   a ⋅cos(b)
    cos(b) - a⋅sin(b) - ───────── + ───────── + ───────── - ───────── - ─────────
                            2           6           24         120         720
    <BLANKLINE>
       7           8           9
      a ⋅sin(b)   a ⋅cos(b)   a ⋅sin(b)    ⎛ 10⎞
    + ───────── + ───────── - ───────── + O⎝a  ⎠
         5040       40320       362880
    >>> g = collect(g.removeO(), [sin(b), cos(b)])
    >>> g
    ⎛   8      6    4    2    ⎞          ⎛     9       7      5    3    ⎞
    ⎜  a      a    a    a     ⎟          ⎜    a       a      a    a     ⎟
    ⎜───── - ─── + ── - ── + 1⎟⋅cos(b) + ⎜- ────── + ──── - ─── + ── - a⎟⋅sin(b)
    ⎝40320   720   24   2     ⎠          ⎝  362880   5040   120   6     ⎠
    >>> sina = sin(a).series(a, 0, 10).removeO()
    >>> sina
       9       7      5    3
      a       a      a    a
    ────── - ──── + ─── - ── + a
    362880   5040   120   6
    >>> # Note that subs will not work directly with sina, because -sina appears
    >>> # in g.  We get around it by substituting -sina.
    >>> h = g.subs(-sina, -sin(a))
    >>> h
    ⎛   8      6    4    2    ⎞
    ⎜  a      a    a    a     ⎟
    ⎜───── - ─── + ── - ── + 1⎟⋅cos(b) - sin(a)⋅sin(b)
    ⎝40320   720   24   2     ⎠
    >>> cosa = cos(a).series(a, 0, 10).removeO()
    >>> cosa
       8      6    4    2
      a      a    a    a
    ───── - ─── + ── - ── + 1
    40320   720   24   2
    >>> h = h.subs(cosa, cos(a))
    >>> h
    -sin(a)⋅sin(b) + cos(a)⋅cos(b)
    >>> eq = Eq(f, h)
    >>> eq
    cos(a + b) = -sin(a)⋅sin(b) + cos(a)⋅cos(b)
    >>> eq.lhs.expand(trig=True) == eq.rhs
    True

Not only symbolics: numerical computing
---------------------------------------

.. _solution_numerics_1:

Solution 1
~~~~~~~~~~

::

    >>> from sympy import log, pi

    >>> f = x**(1 - log(log(log(log(1/x)))))

    >>> f.subs(x, pi).evalf()
    0.730290019485505 - 1.30803618190213⋅ⅈ

.. _solution_numerics_2:

Solution 2
~~~~~~~~~~

::

    >>> from sympy import E, pi

    >>> E_million_digits = str(E.evalf(n=1000000))
    >>> pi_million_digits = str(pi.evalf(n=1000000))
    >>> # We can use the .find() method of str.
    >>> # First, let's see what we should offset the value
    >>> # (it's not immediately obvious because of 0 indexing
    >>> # and because of the '.' in the string.
    >>> E_million_digits[:5]
    2.718
    >>> E_million_digits.find('71')
    2
    >>> # So, counting the 2 as the first digit, we see that
    >>> # there is no need for any offset!
    >>> E999999 = E_million_digits.find('999999')
    >>> E999999
    384341
    >>> "%s" % E_million_digits[E999999:E999999+10]
    9999999951
    >>> # We see that there are actually 8 consecutive 9's here
    >>> pi789 = pi_million_digits.find('789')
    >>> pi789
    352
    >>> "%s" % pi_million_digits[pi789:pi789+10]
    7892590360
    >>> # Now, double check the values from Dinosaur Comics.
    >>> E789 = E_million_digits.find('789')
    >>> E789
    2501
    >>> "%s" % E_million_digits[E789:E789+10]
    7893652605
    >>> pi999999 = pi_million_digits.find('999999')
    >>> pi999999
    763
    >>> # Hey, this one was wrong in the webcomic (T-Rex said "762nd")!
    >>> "%s" % pi_million_digits[pi999999:pi999999+10]
    9999998372

.. _solution_numerics_3:

Solution 3
~~~~~~~~~~

::

    >>> from sympy import limit, erf, exp, oo

    >>> limit((erf(x - exp(-exp(x))) - erf(x))*exp(exp(x))*exp(x**2), x, oo)
      -2
    ─────
      ___
    ╲╱ π

Summing roots of polynomials
----------------------------

.. _solution_rootsum_1:

Solution 1
~~~~~~~~~~

* `f = z^5 + z + a` and `g = \frac{1}{z + 1}`::

    >>> from sympy import viete, expand, symmetrize, Eq, pprint, cancel, RootSum

    >>> var('a,z')
    (a, z)
    >>> f = z**5 + z + a
    >>> f
         5
    a + z  + z
    >>> g = 1/(z + 1)
    >>> g
      1
    ─────
    z + 1
    >>> R = var('r:5')
    >>> G = together(sum( [ g.subs(z, r) for r in R] ))
    >>> G
    (r₀ + 1)⋅(r₁ + 1)⋅(r₂ + 1)⋅(r₃ + 1) + (r₀ + 1)⋅(r₁ + 1)⋅(r₂ + 1)⋅(r₄ + 1) + (r
    ──────────────────────────────────────────────────────────────────────────────
                                                                            (r₀ +
    <BLANKLINE>
    ₀ + 1)⋅(r₁ + 1)⋅(r₃ + 1)⋅(r₄ + 1) + (r₀ + 1)⋅(r₂ + 1)⋅(r₃ + 1)⋅(r₄ + 1) + (r₁
    ──────────────────────────────────────────────────────────────────────────────
    1)⋅(r₁ + 1)⋅(r₂ + 1)⋅(r₃ + 1)⋅(r₄ + 1)
    <BLANKLINE>
    + 1)⋅(r₂ + 1)⋅(r₃ + 1)⋅(r₄ + 1)
    ───────────────────────────────
    <BLANKLINE>
    >>> V = viete(f, R, z)
    >>> for lhs, rhs in V:
    ...     pprint(Eq(lhs, rhs), use_unicode=True)
    ...
    r₀ + r₁ + r₂ + r₃ + r₄ = 0
    r₀⋅r₁ + r₀⋅r₂ + r₀⋅r₃ + r₀⋅r₄ + r₁⋅r₂ + r₁⋅r₃ + r₁⋅r₄ + r₂⋅r₃ + r₂⋅r₄ + r₃⋅r₄
    = 0
    r₀⋅r₁⋅r₂ + r₀⋅r₁⋅r₃ + r₀⋅r₁⋅r₄ + r₀⋅r₂⋅r₃ + r₀⋅r₂⋅r₄ + r₀⋅r₃⋅r₄ + r₁⋅r₂⋅r₃ + r
    ₁⋅r₂⋅r₄ + r₁⋅r₃⋅r₄ + r₂⋅r₃⋅r₄ = 0
    r₀⋅r₁⋅r₂⋅r₃ + r₀⋅r₁⋅r₂⋅r₄ + r₀⋅r₁⋅r₃⋅r₄ + r₀⋅r₂⋅r₃⋅r₄ + r₁⋅r₂⋅r₃⋅r₄ = 1
    r₀⋅r₁⋅r₂⋅r₃⋅r₄ = -a
    >>> p = expand(numer(G))
    >>> q = expand(denom(G))
    >>> (P, Q), mapping = symmetrize((p, q), R, formal=True)
    >>> P
    (4⋅s₁ + 3⋅s₂ + 2⋅s₃ + s₄ + 5, 0)
    >>> Q
    (s₁ + s₂ + s₃ + s₄ + s₅ + 1, 0)
    >>> for s, poly in mapping:
    ...     pprint(Eq(s, poly), use_unicode=True)
    ...
    s₁ = r₀ + r₁ + r₂ + r₃ + r₄
    s₂ = r₀⋅r₁ + r₀⋅r₂ + r₀⋅r₃ + r₀⋅r₄ + r₁⋅r₂ + r₁⋅r₃ + r₁⋅r₄ + r₂⋅r₃ + r₂⋅r₄ + r
    ₃⋅r₄
    s₃ = r₀⋅r₁⋅r₂ + r₀⋅r₁⋅r₃ + r₀⋅r₁⋅r₄ + r₀⋅r₂⋅r₃ + r₀⋅r₂⋅r₄ + r₀⋅r₃⋅r₄ + r₁⋅r₂⋅r
    ₃ + r₁⋅r₂⋅r₄ + r₁⋅r₃⋅r₄ + r₂⋅r₃⋅r₄
    s₄ = r₀⋅r₁⋅r₂⋅r₃ + r₀⋅r₁⋅r₂⋅r₄ + r₀⋅r₁⋅r₃⋅r₄ + r₀⋅r₂⋅r₃⋅r₄ + r₁⋅r₂⋅r₃⋅r₄
    s₅ = r₀⋅r₁⋅r₂⋅r₃⋅r₄
    >>> subslist = [ (s, c) for (s, _), (_, c) in zip(mapping, V) ]
    >>> subslist
    [(s₁, 0), (s₂, 0), (s₃, 0), (s₄, 1), (s₅, -a)]
    >>> P[0].subs(subslist)/Q[0].subs(subslist)
      6
    ──────
    -a + 2
    >>> # Note, we must give a variable to RootSum, as f has more than one Symbol
    >>> cancel(P[0].subs(subslist)/Q[0].subs(subslist) - RootSum(f, Lambda(z, g), z)) == 0
    True


* `f = z^5 + z + a` and `g = \frac{1}{z + b}`::

    >>> var('a b')
    (a, b)
    >>> f = z**5 + z + a
    >>> f
         5
    a + z  + z
    >>> g = 1/(z + b)
    >>> g
      1
    ─────
    b + z
    >>> R = var('r:5')
    >>> G = together(sum( [ g.subs(z, r) for r in R] ))
    >>> G
    (b + r₀)⋅(b + r₁)⋅(b + r₂)⋅(b + r₃) + (b + r₀)⋅(b + r₁)⋅(b + r₂)⋅(b + r₄) + (b
    ──────────────────────────────────────────────────────────────────────────────
                                                                            (b + r
    <BLANKLINE>
     + r₀)⋅(b + r₁)⋅(b + r₃)⋅(b + r₄) + (b + r₀)⋅(b + r₂)⋅(b + r₃)⋅(b + r₄) + (b +
    ──────────────────────────────────────────────────────────────────────────────
    ₀)⋅(b + r₁)⋅(b + r₂)⋅(b + r₃)⋅(b + r₄)
    <BLANKLINE>
     r₁)⋅(b + r₂)⋅(b + r₃)⋅(b + r₄)
    ───────────────────────────────
    <BLANKLINE>
    >>> V = viete(f, R, z)
    >>> for lhs, rhs in V:
    ...     pprint(Eq(lhs, rhs), use_unicode=True)
    ...
    r₀ + r₁ + r₂ + r₃ + r₄ = 0
    r₀⋅r₁ + r₀⋅r₂ + r₀⋅r₃ + r₀⋅r₄ + r₁⋅r₂ + r₁⋅r₃ + r₁⋅r₄ + r₂⋅r₃ + r₂⋅r₄ + r₃⋅r₄
    = 0
    r₀⋅r₁⋅r₂ + r₀⋅r₁⋅r₃ + r₀⋅r₁⋅r₄ + r₀⋅r₂⋅r₃ + r₀⋅r₂⋅r₄ + r₀⋅r₃⋅r₄ + r₁⋅r₂⋅r₃ + r
    ₁⋅r₂⋅r₄ + r₁⋅r₃⋅r₄ + r₂⋅r₃⋅r₄ = 0
    r₀⋅r₁⋅r₂⋅r₃ + r₀⋅r₁⋅r₂⋅r₄ + r₀⋅r₁⋅r₃⋅r₄ + r₀⋅r₂⋅r₃⋅r₄ + r₁⋅r₂⋅r₃⋅r₄ = 1
    r₀⋅r₁⋅r₂⋅r₃⋅r₄ = -a
    >>> p = expand(numer(G))
    >>> q = expand(denom(G))
    >>> (P, Q), mapping = symmetrize((p, q), R, formal=True)
    >>> P
    ⎛   4      3         2                    ⎞
    ⎝5⋅b  + 4⋅b ⋅s₁ + 3⋅b ⋅s₂ + 2⋅b⋅s₃ + s₄, 0⎠
    >>> Q
    ⎛ 5    4       3       2                  ⎞
    ⎝b  + b ⋅s₁ + b ⋅s₂ + b ⋅s₃ + b⋅s₄ + s₅, 0⎠
    >>> for s, poly in mapping:
    ...     pprint(Eq(s, poly), use_unicode=True)
    ...
    s₁ = r₀ + r₁ + r₂ + r₃ + r₄
    s₂ = r₀⋅r₁ + r₀⋅r₂ + r₀⋅r₃ + r₀⋅r₄ + r₁⋅r₂ + r₁⋅r₃ + r₁⋅r₄ + r₂⋅r₃ + r₂⋅r₄ + r
    ₃⋅r₄
    s₃ = r₀⋅r₁⋅r₂ + r₀⋅r₁⋅r₃ + r₀⋅r₁⋅r₄ + r₀⋅r₂⋅r₃ + r₀⋅r₂⋅r₄ + r₀⋅r₃⋅r₄ + r₁⋅r₂⋅r
    ₃ + r₁⋅r₂⋅r₄ + r₁⋅r₃⋅r₄ + r₂⋅r₃⋅r₄
    s₄ = r₀⋅r₁⋅r₂⋅r₃ + r₀⋅r₁⋅r₂⋅r₄ + r₀⋅r₁⋅r₃⋅r₄ + r₀⋅r₂⋅r₃⋅r₄ + r₁⋅r₂⋅r₃⋅r₄
    s₅ = r₀⋅r₁⋅r₂⋅r₃⋅r₄
    >>> subslist = [ (s, c) for (s, _), (_, c) in zip(mapping, V) ]
    >>> subslist
    [(s₁, 0), (s₂, 0), (s₃, 0), (s₄, 1), (s₅, -a)]
    >>> P[0].subs(subslist)/Q[0].subs(subslist)
         4
      5⋅b  + 1
    ───────────
          5
    -a + b  + b
    >>> # Note, we must give a variable to RootSum, as f has more than one Symbol
    >>> cancel(P[0].subs(subslist)/Q[0].subs(subslist) - RootSum(f, Lambda(z, g), z)) == 0
    True

.. _solution_rootsum_2:

Solution 2
~~~~~~~~~~

Vertex `k`--coloring of graphs
------------------------------

.. _solution_colorings_1:

Solution 1
~~~~~~~~~~

.. _solution_colorings_2:

Solution 2
~~~~~~~~~~

Algebraic geometry
------------------

.. _solution_geometric_1:

Solution 1
~~~~~~~~~~

