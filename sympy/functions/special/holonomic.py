"""Holonomic Functions and Differential Operators"""

from __future__ import print_function, division

from sympy import symbols, Symbol, diff, collect, S, Mul
from sympy.polys.polytools import lcm, gcd
from sympy.core.expr import Expr
from sympy.printing import sstr
from sympy.matrices import Matrix
from sympy.core.compatibility import range
from sympy.polys.polytools import DMP


def DiffOperatorAlgebra(base_ring, generator):

    ring = DifferentialOperatorAlgebra(base_ring, generator)
    return (ring, ring.gen)


class DifferentialOperatorAlgebra:

    def __init__(self, base_ring, generator):

        self.base_ring = base_ring
        self.gen = DifferentialOperator([0, 1], self, generator)
        self.gen_str = generator

    def __str__(self):

        string = 'Univariate Differential Operator Algebra in intermediate '\
            + self.gen_str + ' over the base ring ' + \
            (self.base_ring).__str__()

        return string

    __repr__ = __str__


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

    _op_priority = 30

    def __init__(self, list_of_poly, parent_ring, generator=None, var=None):

        if generator is None:
            self.gen_symbol = symbols('Dx', commutative=False)
        else:
            if isinstance(generator, str):
                self.gen_symbol = symbols(generator, commutative=False)
            elif isinstance(generator, Symbol):
                self.gen_symbol = generator

        self.parent_ring = parent_ring

        if isinstance(list_of_poly, list):
            for i, j in enumerate(list_of_poly):
                if isinstance(j, int):
                    list_of_poly[i] = (
                        (self.parent_ring).base_ring).from_sympy(S(j))
                elif not isinstance(j, DMP):
                    list_of_poly[i] = (
                        (self.parent_ring).base_ring).from_sympy(j)

            self.listofpoly = list_of_poly

    def __mul__(self, other):
        """
        Multiplies two DifferentialOperator and returns another
        DifferentialOperator isinstance using the commutation rule
        Dx*a = a*Dx + a'
        """

        ring_0 = self.parent_ring.base_ring.from_sympy(S(0))
        ring_1 = self.parent_ring.base_ring.from_sympy(S(1))

        if isinstance(other, DifferentialOperator):
            if other.listofpoly == [ring_0, ring_1]:
                sol = []
                sol.append(ring_0)
                for i in self.listofpoly:
                    sol.append(i)

                return DifferentialOperator(sol, self.parent_ring, self.gen_symbol)

        gen = DifferentialOperator([0, 1], self.parent_ring, self.gen_symbol)
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

            if not isinstance(other, DMP):
                other = (self.parent_ring.base_ring).from_sympy(other)

            sol = []
            for j in self.listofpoly:
                sol.append(other * j)

            return DifferentialOperator(sol, self.parent_ring, self.gen_symbol)

    def __add__(self, other):

        if isinstance(other, DifferentialOperator):

            list_self = self.listofpoly
            list_other = other.listofpoly

            if min(len(list_self), len(list_other)) is len(list_self):
                minimum = 0
            else:
                minimum = 1

            if minimum is 0:
                sol = [
                    a + b for a, b in zip(list_self, list_other)] + list_other[len(list_self):]
                return DifferentialOperator(sol, self.parent_ring, self.gen_symbol)

            else:
                sol = [
                    a + b for a, b in zip(list_self, list_other)] + list_self[len(list_other):]
                return DifferentialOperator(sol, self.parent_ring, self.gen_symbol)

        else:
            if isinstance(other, int):
                other = S(other)
            list_self = self.listofpoly
            if not isinstance(other, DMP):
                list_other = [((self.parent_ring).base_ring).from_sympy(other)]
            else:
                list_other = [other]
            sol = []
            sol.append(list_self[0] + list_other[0])
            sol += list_self[1:]

            return DifferentialOperator(sol, self.parent_ring, self.gen_symbol)

    __radd__ = __add__

    def __sub__(self, other):

        return self + (-1) * other

    def __rsub__(self, other):
        return (-1) * self + other

    def __pow__(self, n):

        ring_0 = self.parent_ring.base_ring.from_sympy(S(0))
        ring_1 = self.parent_ring.base_ring.from_sympy(S(1))

        if n == 1:
            return self
        if n == 0:
            return DifferentialOperator([1], self.parent_ring, self.gen_symbol)

        if self.listofpoly == [ring_0, ring_1]:
            sol = []
            for i in range(0, n):
                sol.append(self.parent_ring.base_ring.from_sympy(S(0)))
            sol.append(self.parent_ring.base_ring.from_sympy(S(1)))

            return DifferentialOperator(sol, self.parent_ring, self.gen_symbol)

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

    def _normalize(ann):
        """
        Normalize a given annihilator
        """

        if isinstance(ann, DifferentialOperator):
            annihilator = ann
        else:
            annihilator = DifferentialOperator(ann)

        x = annihilator.var
        listofpoly = annihilator.listofpoly
        list_of_coeff = []
        y = symbols('y')

        for i in listofpoly:
            list_of_coeff.append(i.subs(x, y))

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
            list_of_coeff[i] = list_of_coeff[i].subs(y, x)

        return _list_to_ann(list_of_coeff)

    def diff(self):

        listofpoly = self.listofpoly
        sol = []
        for i in listofpoly:
            sol.append(i.diff())
        return DifferentialOperator(sol, self.parent_ring, self.gen_symbol)


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
    Holonomic((1) + (x)Dx, x)


    >>> p = DifferentialOperator(Dx-1)
    >>> q = DifferentialOperator(Dx**2+1)
    >>> HoloFunc(p,x) + HoloFunc(q,x)
    Holonomic((-1) + (1)Dx + (-1)Dx**2 + (1)Dx**3, x)


    >>> p = DifferentialOperator(Dx+1)
    >>> q = DifferentialOperator(Dx**2-1)
    >>> HoloFunc(p,x) + HoloFunc(q,x)
    Holonomic((-1) + (1)Dx**2, x)


    >>> p = DifferentialOperator(x*Dx+1)
    >>> q = DifferentialOperator(Dx+5)
    >>> HoloFunc(p,x) + HoloFunc(q,x)
    Holonomic((10 - 25*x) + (2 - 25*x**2)Dx + (x - 5*x**2)Dx**2, x)


    For details see ore_algebra package in Sage
    """

    def __init__(self, ann, x, *args):
        if len(args) is 2:
            self.cond = args[0]
            self.cond_point = args[1]
        elif len(args) is 1:
            self.cond = args[0]
            self.cond_point = 0
        if isinstance(ann, DifferentialOperator):
            self.ann = ann  # ann is an instance of DifferentialOperator
        else:
            self.ann = DifferentialOperator(ann)
        self.var = x

    def __repr__(self):

        return 'Holonomic(%s, %s)' % ((self.ann).__repr__(), sstr((self.ann).var))

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
        sol = _list_to_ann(sol)

        sol = DifferentialOperator(sol) * (self.ann)

        return HoloFunc(sol._normalize(), (self.ann).var)

    def integrate(self):

        D = DifferentialOperator(self.ann.gen)
        return HoloFunc(self.ann * D, self.var)
