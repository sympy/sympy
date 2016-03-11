"""Holonomic Functions"""

from __future__ import print_function, division
from sympy import symbols, poly, diff, Symbol, sstr, pretty_print, printing, Function, dsolve

class DifferentialOperator(object):
    """
    Represents a Differential Operator in which the generator is at RIGHT, supports
    addition and multiplication.

    Examples
    ========

    >>> from sympy.functions.special.holonomic import DifferentialOperator
    >>> Dx = symbols('Dx')
    >>> DifferentialOperator(x)*DifferentialOperator(Dx+1)
    Dx*x + x
    >>> DifferentialOperator(Dx-1)*DifferentialOperator(x*Dx+1)
    Dx**2*x - Dx*x + 2*Dx - 1

    """

    def __init__(self, equation):
        Dx = Symbol('Dx')
        self.eq = equation

    def __mul__(self, right):
        x, Dx = symbols('x, Dx')
        eq = poly(self.eq)
        if not Dx in poly(eq).gens:
            deg = 0
        else:
            deg = poly(eq).degree('Dx')
        poly_comm = right.eq
        init_coeff = poly(eq, Dx).nth(0)
        sol = init_coeff*poly_comm
        def noncomm_mul(b):
            return b*Dx + b.diff('x')
        for i in xrange(1, deg+1):
            poly_comm = noncomm_mul(poly_comm)
            sol = sol + poly(eq, Dx).nth(i)*poly_comm
        return DifferentialOperator(sol.expand())

    def __add__(self, right):
        return DifferentialOperator(self.eq + right.eq)

    def __repr__(self):
        return sstr(self.eq)

    def __str__(self):
        return sstr(self.eq)

class HoloFunc(object):
    """Represents a Holonomic Function,
    first parameter is the annihilator of the holonomic function,
    second is the variable for the function, initial conditions are optional.

    Addition is implemented and works for some cases.

    Examples
    ========
    >>> Dx = symbols('Dx')
    >>> HoloFunc(x*Dx+1,x)
    holonomic(Dx*x + 1)

    >>> HoloFunc(Dx-1,x) + HoloFunc(Dx**2+1,x)
    holonomic(-Dx**3 + Dx**2 - Dx + 1)

    >>> HoloFunc(Dx+1,x) + HoloFunc(Dx**2-1,x)
    holonomic(-Dx**2 + 1)

    >>> HoloFunc(x*Dx-1,x) + HoloFunc(x*Dx+1,x)
    ### Output should be holonomic(x**2*Dx**2 + x*Dx - 1), normalization is not implemented yet
    holonomic(-Dx**2*x - Dx + 1/x)

    For details see ore_algebra package in Sage
    """

    def __init__(self, annihilator, x, *args):
        if len(args) is 2:
            self.cond = args[0]
            self.cond_point = args[1]
        elif len(args) is 1:
            self.cond = args[0]
            self.cond_point = 0
        self.annihilator = annihilator
        Dx = symbols('Dx')
        self.generator = Dx
        self.var = x

        
    def __repr__(self):
        return 'holonomic(%s)' % sstr(self.annihilator)

    def __add__(self, other):

        deg1 = poly(self.annihilator, self.generator).degree()
        deg2 = poly(other.annihilator, other.generator).degree()
        dim = max(deg1, deg2)
        rowsself = [self.annihilator]; rowsother = [other.annihilator]
        for i in range(dim - deg1):
            derivative1 = (DifferentialOperator(self.generator)*DifferentialOperator(rowsself[-1])).eq
            rowsself.append(derivative1)
        for i in range(dim - deg2):
            derivative2 = (DifferentialOperator(other.generator)*DifferentialOperator(rowsother[-1])).eq
            rowsother.append(derivative2)
        row = rowsself + rowsother
        r = [[poly(expr, self.generator).nth(i) for i in xrange(dim+1)] for expr in (row)]
        homogeneous = [[0 for q in xrange(dim+1)]]
        from sympy.matrices import Matrix
        solcomp = (Matrix(r).transpose()).gauss_jordan_solve(Matrix(homogeneous).transpose())
        sol = (Matrix(r).transpose()).gauss_jordan_solve(Matrix(homogeneous).transpose())[0]
        while sol.is_zero:
            dim += 1
            derivative1 = (DifferentialOperator(self.generator)*DifferentialOperator(rowsself[-1])).eq
            rowsself.append(derivative1)
            derivative2 = (DifferentialOperator(other.generator)*DifferentialOperator(rowsother[-1])).eq
            rowsother.append(derivative2)
            row = rowsself + rowsother
            r = [[poly(expr, self.generator).nth(i) for i in xrange(dim+1)] for expr in (row)]
            homogeneous = [[0 for q in xrange(dim+1)]]
            solcomp = (Matrix(r).transpose()).gauss_jordan_solve(Matrix(homogeneous).transpose())
            sol = (Matrix(r).transpose()).gauss_jordan_solve(Matrix(homogeneous).transpose())[0]
        if sol.is_symbolic():
            sol = solcomp[0]/solcomp[1]
        sol = sol[:dim + 1 - deg1]
        solution = 0
        for i, j in enumerate(sol):
            solution += ((self.generator)**i)*j
        return(HoloFunc((DifferentialOperator(solution)*DifferentialOperator(self.annihilator)), self.var))

        

