"""Holonomic Functions and Differential Operators"""

from __future__ import print_function, division

from sympy import symbols, Symbol, diff, S, Mul
from sympy.polys.polytools import lcm, gcd
from sympy.core.expr import Expr
from sympy.printing import sstr
from sympy.matrices import Matrix
from sympy.core.compatibility import range
from sympy.polys.polytools import DMP


def DiffOperatorAlgebra(base, generator):
    """ A function to create a Differential Operator Algebra.
    The first arguments need to be the base polynomial ring for the algebra
    and the second argument must be a generator.
    Returns the class representing the operator algebra and
    the standard operator for differentiation i.e. the `Dx` operator.

    Examples
    =======
    from sympy.polys.domains import ZZ
    R, Dx = DiffOperatorAlgebra
    """

    ring = DifferentialOperatorAlgebra(base, generator)
    return (ring, ring.diff_operator)


class DifferentialOperatorAlgebra:
    """ The class representing Differential Operator Algebra.
    Defined by the base polynomial ring and the generator.

    The attribute diff_operator is the operator `Dx` which when acted
    on a function does differentiation.
    """

    def __init__(self, base, generator):

        self.base = base

        self.diff_operator = DifferentialOperator(
            [base.zero, base.one], self, generator)

        if generator is None:
            self.gen_symbol = symbols('Dx', commutative=False)
        else:
            if isinstance(generator, str):
                self.gen_symbol = symbols(generator, commutative=False)
            elif isinstance(generator, Symbol):
                self.gen_symbol = generator



    def __str__(self):

        string = 'Univariate Differential Operator Algebra in intermediate '\
            + sstr(self.gen_symbol) + ' over the base ring ' + \
            (self.base).__str__()

        return string

    __repr__ = __str__


class DifferentialOperator(Expr):
    """
    The operator can be created by providing a list of polynomials for each power of Dx
    and the parent ring which must be an instance of DifferentialOperatorAlgebra.

    An Operator can also be created easily using the operator `Dx`. See examples below.

    Examples
    ========

    >>> from sympy.functions.special.holonomic import DifferentialOperator, DiffOperatorAlgebra
    >>> from sympy.polys.domains import ZZ, QQ
    >>> from sympy import symbols
    >>> x = symbols('x')
    >>> R, Dx = DiffOperatorAlgebra(ZZ.old_poly_ring(x),'Dx')

    >>> DifferentialOperator([0,1,x**2], R)
    (0) + (1)Dx + (x**2)Dx**2

    One can use the `Dx` and represent differential operators more conveniently this way also.
    >>> (x*Dx*x + 1 - Dx**2)**2
    (2*x**2 + 2*x + 1) + (4*x**3 + 2*x**2 - 4)Dx + (x**4 - 6*x - 2)Dx**2 + (-2*x**2)Dx**3 + (1)Dx**4

    """

    _op_priority = 30

    def __init__(self, list_of_poly, parent, generator=None, variable=None):

        self.parent = parent
        if variable is None:
            self.variable = symbols('x')
        else:
            self.variable = variable

        if isinstance(list_of_poly, list):
            for i, j in enumerate(list_of_poly):
                if isinstance(j, int):
                    list_of_poly[i] = (
                        (self.parent).base).from_sympy(S(j))
                elif not isinstance(j, self.parent.base.dtype):
                    list_of_poly[i] = (
                        (self.parent).base).from_sympy(j)

            self.listofpoly = list_of_poly

    def args(self):

        return self.gen_symbol

    def __mul__(self, other):
        """
        Multiplies two DifferentialOperator and returns another
        DifferentialOperator isinstance using the commutation rule
        Dx*a = a*Dx + a'
        """

        if isinstance(other, DifferentialOperator):
            if other.listofpoly == other.parent.diff_operator.listofpoly:
                sol = []
                sol.append(other.parent.base.zero)
                for i in self.listofpoly:
                    sol.append(i)

                return DifferentialOperator(sol, self.parent)

        gen = self.parent.diff_operator
        listofpoly = self.listofpoly
        sol = (listofpoly[0] * other)

        def _diff_n_times(b):
            expr = b * gen + b.diff()
            return expr

        for i in range(1, len(listofpoly)):
            other = _diff_n_times(other)
            sol += (listofpoly[i] * other)

        return sol

    def __rmul__(self, other):

        if not isinstance(other, DifferentialOperator):

            if isinstance(other, int):
                other = S(other)

            if not isinstance(other, self.parent.base.dtype):
                other = (self.parent.base).from_sympy(other)

            sol = []
            for j in self.listofpoly:
                sol.append(other * j)

            return DifferentialOperator(sol, self.parent)

    def __add__(self, other):

        if isinstance(other, DifferentialOperator):

            list_self = self.listofpoly
            list_other = other.listofpoly

            if len(list_self) <= len(list_other):
                sol = [
                    a + b for a, b in zip(list_self, list_other)] + list_other[len(list_self):]
                return DifferentialOperator(sol, self.parent)

            else:
                sol = [
                    a + b for a, b in zip(list_self, list_other)] + list_self[len(list_other):]
                return DifferentialOperator(sol, self.parent)

        else:

            if isinstance(other, int):
                other = S(other)
            list_self = self.listofpoly
            if not isinstance(other, self.parent.base.dtype):
                list_other = [((self.parent).base).from_sympy(other)]
            else:
                list_other = [other]
            sol = []
            sol.append(list_self[0] + list_other[0])
            sol += list_self[1:]

            return DifferentialOperator(sol, self.parent)

    __radd__ = __add__

    def __sub__(self, other):

        return self + (-1) * other

    def __rsub__(self, other):
        return (-1) * self + other

    def __pow__(self, n):

        if n == 1:
            return self
        if n == 0:
            return DifferentialOperator([1], self.parent)

        if self.listofpoly == self.parent.diff_operator.listofpoly:
            sol = []
            for i in range(0, n):
                sol.append(self.parent.base.zero)
            sol.append(self.parent.base.one)

            return DifferentialOperator(sol, self.parent)

        else:
            if n % 2 == 1:
                powreduce = self**(n - 1)
                return powreduce * self
            elif n % 2 == 0:
                powreduce = self**(n / 2)
                return powreduce * powreduce

    def __str__(self):

        listofpoly = self.listofpoly
        print_str = ''

        for i, j in enumerate(listofpoly):
            if j == 0:
                continue

            if i == 0:
                print_str += '(' + sstr(j) + ')'
                continue

            if print_str:
                print_str += ' + '

            if i == 1:
                print_str += '(' + sstr(j) + ')Dx'
                continue

            print_str += '(' + sstr(j) + ')' + 'Dx**' + sstr(i)

        return print_str

    __repr__ = __str__

    def diff(self):

        listofpoly = self.listofpoly
        sol = []
        for i in listofpoly:
            sol.append(i.diff())
        return DifferentialOperator(sol, self.parent)

def _normalize(list_of_coeff, x):
    """
    Normalize a given annihilator
    """

    lcm_denom = S(1)

    for i in list_of_coeff:
        lcm_denom = lcm(lcm_denom, i.as_numer_denom()[1])

    if isinstance(list_of_coeff[-1], Mul) and (list_of_coeff[-1].args[0]).is_negative:
        lcm_denom = -lcm_denom
    elif list_of_coeff[-1].is_negative:
        lcm_denom = -lcm_denom

    for i, j in enumerate(list_of_coeff):
        list_of_coeff[i] = (j * lcm_denom).simplify()

    gcd_numer = list_of_coeff[-1]

    for i in list_of_coeff:
        gcd_numer = gcd(gcd_numer, i.as_numer_denom()[0])

    for i, j in enumerate(list_of_coeff):
        list_of_coeff[i] = (j / gcd_numer).simplify()

    return list_of_coeff

class HoloFunc(object):
    """Represents a Holonomic Function,
    first parameter is the annihilator of the holonomic function, which
    is an instance of DifferentialOperator, second is the variable
    for the function, initial conditions are optional.

    Addition is implemented and works for some cases.

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
        self.var = x

    def __repr__(self):

        return 'Holonomic(%s, %s)' % ((self.annihilator).__repr__(), sstr(self.var))

    def __add__(self, other):

        deg1 = len((self.annihilator).listofpoly) - 1
        deg2 = len((other.annihilator).listofpoly) - 1
        dim = max(deg1, deg2)

        rowsself = [self.annihilator]
        rowsother = [other.annihilator]
        gen = self.annihilator.parent.diff_operator

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
                    temp = self.annihilator.parent.base.to_sympy(
                        expr.listofpoly[i])
                    p.append(temp)
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
                        temp = self.annihilator.parent.base.to_sympy(
                            expr.listofpoly[i])
                        p.append(temp)
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
        sol = _normalize(sol, self.var)
        # construct expression from the coefficients
        sol1 = DifferentialOperator(
            sol, self.annihilator.parent, self.var)

        sol = sol1 * (self.annihilator)

        return HoloFunc(sol, self.var)

    def integrate(self):

        D = DifferentialOperator(self.ann.gen)
        return HoloFunc(self.ann * D, self.var)
