"""Holonomic Functions and Differential Operators"""

from __future__ import print_function, division

from sympy import symbols, poly, diff, Symbol, sstr, pretty_print, printing, Function, Expr, collect, S
from sympy.matrices import Matrix
from sympy.core.compatibility import range


def list_to_ann(alist, gen=None):
    """
    A function to convert a list of coefficient polynomials
    to the Operator
    """

    if gen is None:
        gen = symbols('Dx', commutative=False)

    expr = 0
    for i, j in enumerate(alist):
        expr += j * gen**i

    return expr


class DifferentialOperator(Expr):
    """
    Represents a Differential Operator whose base ring is the polynomial ring, and supports
    addition and multiplication.

    Examples
    ========

    >>> from sympy import symbols
    >>> from sympy.functions.special.holonomic import DifferentialOperator, HoloFunc
    >>> Dx, x = symbols('Dx, x', commutative=False)
    >>> DifferentialOperator(Dx*x)*DifferentialOperator(Dx+1)
    DifferentialOperator(1 + Dx + x*Dx + x*Dx**2)
    >>> DifferentialOperator(Dx-1)*DifferentialOperator(Dx*x**2+1)
    DifferentialOperator(1 + Dx - 2*x + 4*x*Dx - x**2*Dx + x**2*Dx**2)

    """

    def __init__(self, expr, generator=None, var=None):
        if generator is None:
            generator = symbols('Dx', commutative=False)
        if var is None:
            var = symbols('x', commutative=False)
        self.gen = generator
        self.var = var
        self.expr = expr
        self.expr = self.shift_right()
        self.listofpoly = self.list_of_poly()

    def __mul__(self, other):
        """
        Multiplies two DifferentialOperator and returns another
        DifferentialOperator isinstance using the commutation rule
        Dx*a = a*Dx + a'
        """

        listofpoly = self.listofpoly

        if isinstance(other, DifferentialOperator):
            other = other.expr
        sol = listofpoly[0] * other
        gen = self.gen

        def diffntimes(b):
            expr = b * gen + b.diff(self.var)
            return expr

        for i in range(1, len(listofpoly)):
            other = diffntimes(other)
            sol += listofpoly[i] * other

        return DifferentialOperator(sol.expand())

    def __add__(self, other):
        expr = self.expr + other.expr
        return DifferentialOperator(expr.expand())

    def shift_right(self):
        """
        Substitutes Dx*x --> x*Dx + 1 until all the generators
        comes to right.
        """

        expr = self.expr
        generator = self.gen
        var = self.var
        toright = None

        while not toright is expr:
            toright = expr
            expr = expr.subs(generator * var, var * generator + 1).expand()

        return expr

    def __str__(self):
        return sstr(self.expr)

    __repr__ = __str__

    def list_of_poly(self):
        """
        Converts the noncommutative annihilator expression
        into a list of polynomial for each power of Dx
        """

        gen = self.gen
        dict_coeff = collect(self.expr, gen, evaluate=False)
        listofpoly = []
        r = S(1)

        while not len(dict_coeff) == 0:
            if r in dict_coeff.keys():
                listofpoly.append(dict_coeff[r])
                dict_coeff.pop(r)
            else:
                listofpoly.append(S(0))
            r = r * gen
        return listofpoly


class HoloFunc(object):
    """Represents a Holonomic Function,
    first parameter is the annihilator of the holonomic function, which
    is an instance of DifferentialOperator, second is the variable
    for the function, initial conditions are optional.

    Addition is implemented and works for some cases.

    Examples
    ========
    >>> from sympy import symbols
    >>> from sympy.functions.special.holonomic import DifferentialOperator, HoloFunc
    >>> Dx, x = symbols('Dx, x', commutative=False)
    >>> p = DifferentialOperator(x*Dx+1)
    >>> HoloFunc(p,x)
    holonomic(1 + x*Dx, x)


    >>> p = DifferentialOperator(Dx-1)
    >>> q = DifferentialOperator(Dx**2+1)
    >>> HoloFunc(p,x) + HoloFunc(q,x)
    holonomic(1 - Dx + Dx**2 - Dx**3, x)


    >>> p = DifferentialOperator(Dx+1)
    >>> q = DifferentialOperator(Dx**2-1)
    >>> HoloFunc(p,x) + HoloFunc(q,x)
    holonomic(1 - Dx**2, x)


    # Output differs by a factor
    #normalization is not implemented yet, but its still correct

    >>> p = DifferentialOperator(x*Dx+1)
    >>> q = DifferentialOperator(x*Dx+5)
    >>> HoloFunc(p,x) + HoloFunc(q,x)
    holonomic(-7*Dx - 5*x**(-1) - x*Dx**2, x)


    For details see ore_algebra package in Sage
    """

    def __init__(self, ann, x, *args):
        if len(args) is 2:
            self.cond = args[0]
            self.cond_point = args[1]
        elif len(args) is 1:
            self.cond = args[0]
            self.cond_point = 0
        self.ann = ann  # ann is an instance of DifferentialOperator

    def __repr__(self):

        return 'holonomic(%s, %s)' % ((self.ann).__repr__(), sstr((self.ann).var))

    def __add__(self, other):

        deg1 = len((self.ann).listofpoly) - 1
        deg2 = len((other.ann).listofpoly) - 1
        dim = max(deg1, deg2)

        rowsself = [self.ann]
        rowsother = [other.ann]

        gen = DifferentialOperator((self.ann).gen)

        # constructing annihilators up to order dim
        for i in range(dim - deg1):
            diff1 = (gen * rowsself[-1])
            rowsself.append(diff1)

        for i in range(dim - deg2):
            diff2 = (gen * rowsother[-1])
            rowsother.append(diff2)

        row = rowsself + rowsother

        # constructing the matrix of the ansatz
        r = []

        for expr in row:
            p = []
            for i in range(dim + 1):
                if i >= len(expr.listofpoly):
                    p.append(0)
                else:
                    p.append(expr.listofpoly[i])
            r.append(p)
        r = Matrix(r).transpose()

        homosys = [[0 for q in range(dim + 1)]]
        homosys = Matrix(homosys).transpose()

        # solving the linear system using gauss jordan solver
        solcomp = r.gauss_jordan_solve(homosys)
        sol = (r).gauss_jordan_solve(homosys)[0]

        # if a solution is not obtained then increasing the order by 1 in each
        # iteration
        while sol.is_zero:
            dim += 1

            diff1 = (gen * rowsself[-1])
            rowsself.append(diff1)

            diff2 = (gen * rowsother[-1])
            rowsother.append(diff2)

            row = rowsself + rowsother
            r = []

            for expr in row:
                p = []
                for i in range(dim + 1):
                    if i >= len(expr.listofpoly):
                        p.append(0)
                    else:
                        p.append(expr.listofpoly[i])
                r.append(p)

            r = Matrix(r).transpose()

            homosys = [[0 for q in range(dim + 1)]]
            homosys = Matrix(homosys).transpose()

            solcomp = r.gauss_jordan_solve(homosys)
            sol = r.gauss_jordan_solve(homosys)[0]

        # removing the symbol if any from the solution
        if sol.is_symbolic():
            sol = solcomp[0] / solcomp[1]

        # taking only the coefficients needed to multiply with `self`
        # can be also be done the other way by taking R.H.S and multiply with
        # `other`
        sol = sol[:dim + 1 - deg1]

        # construct expression from the coefficients
        sol = list_to_ann(sol)

        return HoloFunc(DifferentialOperator(sol) * (self.ann), (self.ann).var)
