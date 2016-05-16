"""Holonomic Functions and Differential Operators"""

from __future__ import print_function, division

from sympy import symbols, Symbol, diff, S, Mul
from sympy.polys.polytools import lcm, gcd
from sympy.core.expr import Expr
from sympy.printing import sstr
from sympy.matrices import Matrix
from sympy.core.compatibility import range
from sympy.polys.polytools import DMP
from sympy.functions.combinatorial.factorials import binomial
from sympy.core.sympify import sympify


def DiffOperatorAlgebra(base, generator):
    """
    Returns an Algebra of Differential Operators and the operator for
    differentiation i.e. the `Dx` operator.
    The first argument needs to be the base polynomial ring for the algebra
    and the second argument must be a generator which can be either a
    noncommutative Symbol or a string.

    Examples
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy import symbols
    >>> from sympy.holonomic.holonomic import DiffOperatorAlgebra
    >>> x = symbols('x')
    >>> R, Dx = DiffOperatorAlgebra(ZZ.old_poly_ring(x), 'Dx')
    """

    ring = DifferentialOperatorAlgebra(base, generator)
    return (ring, ring.derivative_operator)


class DifferentialOperatorAlgebra(object):
    """
    An Ore Algebra is a set of noncommutative polynomials in the
    intermediate `Dx` and coefficients in a base ring A. It follows the
    commutation rule:
    Dx * a = sigma(a) * Dx + delta(a)

    Where sigma: A --> A is an endomorphism and delta: A --> A is a
    skew-derivation i.e. delta(ab) = delta(a) * b + sigma(a) * delta(b)

    If one takes the sigma as identity map and delta as the standard derivation
    then it becomes the algebra of Differential Operators also called
    a Weyl Algebra i.e. an algebra whose elements are Differential Operators.

    This class represents a Weyl Algebra and serves as the parent ring for
    Differential Operators.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy import symbols
    >>> from sympy.holonomic.holonomic import DiffOperatorAlgebra
    >>> x = symbols('x')
    >>> R, Dx = DiffOperatorAlgebra(ZZ.old_poly_ring(x), 'Dx')
    >>> R
    Univariate Differential Operator Algebra in intermediate Dx over the base ring
    ZZ[x]

    See Also
    ========

    DifferentialOperator
    """

    def __init__(self, base, generator):
        # the base ring for the algebra
        self.base = base
        # the operator representing differentiation i.e. `Dx`
        self.derivative_operator = DifferentialOperator(
            [base.zero, base.one], self)

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


def _add_lists(list1, list2):
    if len(list1) <= len(list2):
        sol = [a + b for a, b in zip(list1, list2)] + list2[len(list1):]
    else:
        sol = [a + b for a, b in zip(list1, list2)] + list1[len(list2):]
    return sol


class DifferentialOperator(object):
    """
    Differential Operators are elements of Weyl Algebra. The Operators
    are defined by a list of polynomials in the base ring and the
    parent ring of the Operator.

    Takes a list of polynomials for each power of Dx and the
    parent ring which must be an instance of DifferentialOperatorAlgebra.

    A Differential Operator can be created easily using
    the operator `Dx`. See examples below.

    Examples
    ========

    >>> from sympy.holonomic.holonomic import DifferentialOperator, DiffOperatorAlgebra
    >>> from sympy.polys.domains import ZZ, QQ
    >>> from sympy import symbols
    >>> x = symbols('x')
    >>> R, Dx = DiffOperatorAlgebra(ZZ.old_poly_ring(x),'Dx')

    >>> DifferentialOperator([0, 1, x**2], R)
    (1)Dx + (x**2)Dx**2

    >>> (x*Dx*x + 1 - Dx**2)**2
    (2*x**2 + 2*x + 1) + (4*x**3 + 2*x**2 - 4)Dx + (x**4 - 6*x - 2)Dx**2 + (-2*x**2)Dx**3 + (1)Dx**4

    See Also
    ========

    DifferentialOperatorAlgebra
    """

    _op_priority = 20

    def __init__(self, list_of_poly, parent):
        # the parent ring for this operator
        # must be an DifferentialOperatorAlgebra object
        self.parent = parent
        # sequence of polynomials in x for each power of Dx
        # represents the operator
        # convert the expressions into ring elements using from_sympy
        if isinstance(list_of_poly, list):
            for i, j in enumerate(list_of_poly):
                if isinstance(j, int):
                    list_of_poly[i] = (
                        (self.parent).base).from_sympy(S(j))
                elif not isinstance(j, self.parent.base.dtype):
                    list_of_poly[i] = (
                        (self.parent).base).from_sympy(j)

            self.listofpoly = list_of_poly
        self.order = len(self.listofpoly) - 1

    def __mul__(self, other):
        """
        Multiplies two DifferentialOperator and returns another
        DifferentialOperator instance using the commutation rule
        Dx*a = a*Dx + a'
        """

        listofself = self.listofpoly

        if not isinstance(other, DifferentialOperator):
            if not isinstance(other, self.parent.base.dtype):
                listofother = [self.parent.base.from_sympy(other)]

            else:
                listofother = [other]
        else:
            listofother = other.listofpoly
        # multiply a polynomial `b` with a list of polynomials

        def _mul_dmp_diffop(b, listofother):

            if isinstance(listofother, list):
                sol = []
                for i in listofother:
                    sol.append(i * b)
                return sol
            else:
                return [b * listofother]

        sol = _mul_dmp_diffop(listofself[0], listofother)

        # compute Dx^i * b
        def _mul_Dxi_b(b):

            sol1 = [self.parent.base.zero]
            sol2 = []

            if isinstance(b, list):
                for i in b:
                    sol1.append(i)
                    sol2.append(i.diff())
            else:
                sol1.append(self.parent.base.from_sympy(b))
                sol2.append(self.parent.base.from_sympy(b).diff())

            return _add_lists(sol1, sol2)

        for i in range(1, len(listofself)):
            # find Dx^i * b in ith iteration
            listofother = _mul_Dxi_b(listofother)
            # solution = solution + listofself[i] * (Dx^i * b)
            sol = _add_lists(sol, _mul_dmp_diffop(listofself[i], listofother))

        return DifferentialOperator(sol, self.parent)

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

            sol = _add_lists(self.listofpoly, other.listofpoly)
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
            return DifferentialOperator([self.parent.base.one], self.parent)
        # if self is `Dx`
        if self.listofpoly == self.parent.derivative_operator.listofpoly:
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
            if j == self.parent.base.zero:
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

    def __eq__(self, other):

        if isinstance(other, DifferentialOperator):
            if self.listofpoly == other.listofpoly and self.parent == other.parent:
                return True
            else:
                return False
        else:
            if self.listofpoly[0] == other:
                for i in listofpoly[1:]:
                    if i is not self.parent.base.zero:
                        return False
                return True
            else:
                return False


def _normalize(list_of_coeff, parent, negative=True):
    """
    Normalize a given annihilator
    """

    num = []
    denom = []
    base = parent.base
    lcm_denom = base.from_sympy(S(1))
    K = base.get_field()

    for i, j in enumerate(list_of_coeff):
        list_of_coeff[i] = K.from_sympy(sympify(j))
        num.append(base(list_of_coeff[i].num))
        denom.append(base(list_of_coeff[i].den))


    for i in denom:
        lcm_denom = i.lcm(lcm_denom)

    if negative is True:
        lcm_denom = -lcm_denom

    for i, j in enumerate(list_of_coeff):
        list_of_coeff[i] = j * K.new(lcm_denom.rep)

    if len(list_of_coeff[-1].den) == 1 and list_of_coeff[-1].den[0] != 1:
        for i, j in enumerate(list_of_coeff):
            list_of_coeff[i] = j * list_of_coeff[-1].den[0]
    gcd_numer = base.from_FractionField(list_of_coeff[-1], K)

    for i in num:
        gcd_numer = i.gcd(gcd_numer)

    for i, j in enumerate(list_of_coeff):
        list_of_coeff[i] = base.from_FractionField(j / K.new(gcd_numer.rep), K)

    return DifferentialOperator(list_of_coeff, parent)


def _derivate_diff_eq(listofpoly):
    """
    Let a differential equation a0(x)y(x) + a1(x)y'(x) + ... = 0
    where a0, a1,... are polynomials or rational functions. The function
    returns b0, b1, b2... such that the differential equation
    b0(x)y(x) + b1(x)y'(x) +... = 0 is formed after differentiating the
    former equation.
    """

    sol = []
    a = len(listofpoly) - 1
    sol.append(listofpoly[0].diff())

    for i, j in enumerate(listofpoly[1:]):
        sol.append(j.diff() + listofpoly[i])

    sol.append(listofpoly[a])
    return sol


class HolonomicFunction(object):
    """
    A Holonomic Function is a solution to a linear homogeneous ordinary
    differential equation with polynomial coefficients. This differential
    equation can also be represented by an annihilator i.e. a Differential
    Operator L such that L.f = 0. For uniqueness of these functions,
    initial conditions can also be provided along with the annihilator.

    Holonomic functions have closure properties and thus forms a ring.
    Given two Holonomic Functions f and g, their sum, product,
    integral and derivative is also a Holonomic Function.

    Examples
    ========

    >>> from sympy.holonomic.holonomic import HolonomicFunction, DiffOperatorAlgebra
    >>> from sympy.polys.domains import ZZ, QQ
    >>> from sympy import symbols
    >>> x = symbols('x')
    >>> R, Dx = DiffOperatorAlgebra(QQ.old_poly_ring(x),'Dx')

    >>> p = HolonomicFunction(Dx - 1, x, 0, [1]) # e^x
    >>> q = HolonomicFunction(Dx**2 + 1, x, 0, [0, 1]) # sin(x)

    >>> p + q # annihilator of e^x + sin(x)
    HolonomicFunction((-1) + (1)Dx + (-1)Dx**2 + (1)Dx**3, x), f(0) = 1 , f'(0) = 2 , f''(0) = 1

    >>> p * q # annihilator of e^x * sin(x)
    HolonomicFunction((2) + (-2)Dx + (1)Dx**2, x), f(0) = 0 , f'(0) = 1 , f''(0) = 2
    """

    def __init__(self, annihilator, x, x0=0, y0=[]):
        """
        Takes the annihilator and variable of the function.
        x0 is the point for which initial conditions are given and
        y0 is a vector of initial values y0 = [f(x0), f'(x0), f''(x0) ...]
        """

        if not isinstance(y0, list):
            self.y0 = [y0]
        else:
            self.y0 = y0

        if len(self.y0) == 0:
            self._have_init_cond = False
        else:
            self._have_init_cond = True

        self.x0 = x0
        # differential operator L such that L.f = 0
        self.annihilator = annihilator
        self.x = x

    def __repr__(self):

        str_sol = 'HolonomicFunction(%s, %s)' % ((self.annihilator).__repr__(), sstr(self.x))
        if not self._have_init_cond:
            return str_sol
        else:
            cond_str = ''
            diff_str = ''
            for i in self.y0:
                cond_str += ', f%s(%d) = %d ' %(diff_str, self.x0, i)
                diff_str += "'"

            sol = str_sol + cond_str
            return sol

    __str__ = __repr__

    def __add__(self, other):

        deg1 = self.annihilator.order
        deg2 = other.annihilator.order
        dim = max(deg1, deg2)

        rowsself = [self.annihilator]
        rowsother = [other.annihilator]
        gen = self.annihilator.parent.derivative_operator

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
                    if isinstance(expr.listofpoly[i], self.annihilator.parent.base.dtype):
                        ring_to_expr = self.annihilator.parent.base.to_sympy(
                            expr.listofpoly[i])
                        p.append(ring_to_expr)
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
                        if isinstance(expr.listofpoly[i], self.annihilator.parent.base.dtype):
                            ring_to_expr = self.annihilator.parent.base.to_sympy(
                                expr.listofpoly[i])
                            p.append(ring_to_expr)
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
        sol1 = _normalize(sol, self.annihilator.parent)
        # annihilator of the solution
        sol = sol1 * (self.annihilator)

        # solving initial conditions
        if self._have_init_cond and other._have_init_cond:
            if self.x0 == other.x0:
                # try to extended the initial conditions
                # using the annihilator
                y0_self = _extend_y0(self, sol.order)
                y0_other = _extend_y0(other, sol.order)
                y0 = [a + b for a, b in zip(y0_self, y0_other)]
                return HolonomicFunction(sol, self.x, self.x0, y0)
            else:
                raise NotImplementedError

        return HolonomicFunction(sol, self.x)

    def integrate(self):
        # just multiply by Dx from right
        D = self.annihilator.parent.derivative_operator
        return HolonomicFunction(self.annihilator * D, self.x)

    def __eq__(self, other):

        if self.annihilator == other.annihilator:
            if self.x == other.x:
                if self._have_init_cond and other._have_init_cond:
                    if self.x0 == other.x0 and self.y0 == other.y0:
                        return True
                    else:
                        return False
                else:
                    return True
            else:
                return False
        else:
            return False

    def __mul__(self, other):

        ann_self = self.annihilator
        ann_other = other.annihilator
        list_self = ann_self.listofpoly
        list_other = ann_other.listofpoly
        a = ann_self.order
        b = ann_other.order

        for i, j in enumerate(list_self):
            if isinstance(j, ann_self.parent.base.dtype):
                list_self[i] = ann_self.parent.base.to_sympy(j)

        for i, j in enumerate(list_other):
            if isinstance(j, ann_other.parent.base.dtype):
                list_other[i] = ann_other.parent.base.to_sympy(j)
        # will be used to reduce the degree
        self_red = [-list_self[i] / list_self[a]
                    for i in range(a)]

        other_red = [-list_other[i] / list_other[b]
                     for i in range(b)]
        # coeff_mull[i][j] is the coefficient of Dx^i(f).Dx^j(g)
        coeff_mul = [[S(0) for i in range(b + 1)]
                     for j in range(a + 1)]
        coeff_mul[0][0] = S(1)
        # making the ansatz
        lin_sys = [[coeff_mul[i][j]
                    for i in range(a) for j in range(b)]]

        homo_sys = [[0 for q in range(a * b)]]
        homo_sys = Matrix(homo_sys).transpose()

        sol = (Matrix(lin_sys).transpose()).gauss_jordan_solve(homo_sys)
        # until a non trivial solution is found
        while sol[0].is_zero:
            # updating the coefficents Dx^i(f).Dx^j(g) for next degree
            for i in range(a - 1, -1, -1):
                for j in range(b - 1, -1, -1):
                    coeff_mul[i][j + 1] += coeff_mul[i][j]
                    coeff_mul[i + 1][j] += coeff_mul[i][j]
                    coeff_mul[i][j] = coeff_mul[i][j].diff()
            # reduce the terms to lower power using annihilators of f, g
            for i in range(a + 1):
                if not coeff_mul[i][b] == S(0):
                    for j in range(b):
                        coeff_mul[i][j] += other_red[j] * \
                            coeff_mul[i][b]
                    coeff_mul[i][b] = S(0)

            # not d2 + 1, as that is already covered in previous loop
            for j in range(b):
                if not coeff_mul[a][j] == 0:
                    for i in range(a):
                        coeff_mul[i][j] += self_red[i] * \
                            coeff_mul[a][j]
                    coeff_mul[a][j] = S(0)

            lin_sys.append([coeff_mul[i][j] for i in range(a)
                            for j in range(b)])

            sol = (Matrix(lin_sys).transpose()).gauss_jordan_solve(homo_sys)

        sol = sol[0] / sol[1]

        sol_ann = _normalize(sol[0:], self.annihilator.parent, negative=False)

        if self._have_init_cond and other._have_init_cond:
            if self.x0 == other.x0:

                # try to find more inital conditions
                y0_self = _extend_y0(self, sol_ann.order)
                y0_other = _extend_y0(other, sol_ann.order)
                # h(x0) = f(x0) * g(x0)
                y0 = [y0_self[0] * y0_other[0]]

                # coefficient of Dx^j(f)*Dx^i(g) in Dx^i(fg)
                for i in range(1, min(len(y0_self), len(y0_other))):
                    coeff = [[0 for i in range(i + 1)] for j in range(i + 1)]
                    for j in range(i + 1):
                        for k in range(i + 1):
                            if j + k == i:
                                coeff[j][k] = binomial(i, j)

                    sol = 0
                    for j in range(i + 1):
                        for k in range(i + 1):
                            sol += coeff[j][k]* y0_self[j] * y0_other[k]

                    y0.append(sol)
                return HolonomicFunction(sol_ann, self.x, self.x0, y0)
            else:
                raise NotImplementedError

        return HolonomicFunction(sol_ann, self.x)


def _extend_y0(Holonomic, n):
    """
    Tries to find more initial conditions by substituting the initial
    value point in the differential equation.
    """

    annihilator = Holonomic.annihilator
    a = annihilator.order
    x = Holonomic.x
    listofpoly = annihilator.listofpoly
    y0 = Holonomic.y0
    for i, j in enumerate(listofpoly):
            if isinstance(j, annihilator.parent.base.dtype):
                listofpoly[i] = annihilator.parent.base.to_sympy(j)

    if len(y0) < a or n <= len(y0):
        return y0
    else:
        list_red = [-listofpoly[i] / listofpoly[a]
                    for i in range(a)]
        y1 = y0
        for i in range(n - a):
            sol = 0
            for a, b in zip(y1, list_red):
                if not b.subs(x, Holonomic.x0).is_finite:
                    return y0
                sol += a * b.subs(x, Holonomic.x0)
            y1.append(sol)
            list_red = _derivate_diff_eq(list_red)

        return y0 + y1[len(y0):]
