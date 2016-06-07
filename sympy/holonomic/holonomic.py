"""Holonomic Functions and Differential Operators"""

from __future__ import print_function, division

from sympy import symbols, Symbol, diff, S, Dummy, Order, rf, meijerint
from sympy.printing import sstr
from .linearsolver import NewMatrix
from .recurrence import HolonomicSequence, RecurrenceOperator, RecurrenceOperators
from sympy.core.compatibility import range
from sympy.functions.combinatorial.factorials import binomial, factorial
from sympy.core.sympify import sympify
from sympy.polys.domains import QQ, ZZ
from sympy.polys.domains.pythonrational import PythonRational
from sympy.simplify.hyperexpand import hyperexpand
from sympy.functions.special.hyper import hyper, meijerg
from sympy.core.numbers import NaN, Infinity, NegativeInfinity
from sympy.matrices import Matrix
from sympy.polys.polyclasses import DMF
from sympy.polys.polyroots import roots
from sympy.functions.elementary.exponential import exp_polar, exp


def DifferentialOperators(base, generator):
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
    >>> from sympy.holonomic.holonomic import DifferentialOperators
    >>> x = symbols('x')
    >>> R, Dx = DifferentialOperators(ZZ.old_poly_ring(x), 'Dx')
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
    >>> from sympy.holonomic.holonomic import DifferentialOperators
    >>> x = symbols('x')
    >>> R, Dx = DifferentialOperators(ZZ.old_poly_ring(x), 'Dx')
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

    def __eq__(self, other):
        if self.base == other.base and self.gen_symbol == other.gen_symbol:
            return True
        else:
            return False


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

    >>> from sympy.holonomic.holonomic import DifferentialOperator, DifferentialOperators
    >>> from sympy.polys.domains import ZZ, QQ
    >>> from sympy import symbols
    >>> x = symbols('x')
    >>> R, Dx = DifferentialOperators(ZZ.old_poly_ring(x),'Dx')

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
                    list_of_poly[i] = self.parent.base.from_sympy(S(j))
                elif not isinstance(j, self.parent.base.dtype):
                    list_of_poly[i] = self.parent.base.from_sympy(j)

            self.listofpoly = list_of_poly
        # highest power of `Dx`
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
                listofother = [self.parent.base.from_sympy(sympify(other))]

            else:
                listofother = [other]
        else:
            listofother = other.listofpoly

        # multiplies a polynomial `b` with a list of polynomials
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

        # the general case
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


def _normalize(list_of, parent, negative=True):
    """
    Normalize a given annihilator
    """

    num = []
    denom = []
    base = parent.base
    K = base.get_field()
    R = ZZ.old_poly_ring(base.gens[0])
    lcm_denom = R.from_sympy(S(1))
    list_of_coeff = []

    # convert polynomials to the elements of associated
    # fraction field
    for i, j in enumerate(list_of):
        if not isinstance(j, K.dtype):
            list_of_coeff.append(K.from_sympy(sympify(j)))
        else:
            list_of_coeff.append(j)

        # corresponding numerators of the sequence of polynomials
        num.append(base(list_of_coeff[i].num))

        # corresponding denominators
        den = list_of_coeff[i].den
        if isinstance(den[0], PythonRational):
            for i, j in enumerate(den):
                den[i] = j.p

        denom.append(R(den))

    # lcm of denominators in the coefficients
    for i in denom:
        lcm_denom = i.lcm(lcm_denom)

    if negative is True:
        lcm_denom = -lcm_denom

    lcm_denom = K.new(lcm_denom.rep)

    # multiply the coefficients with lcm
    for i, j in enumerate(list_of_coeff):
        list_of_coeff[i] = j * lcm_denom

    gcd_numer = base.from_FractionField(list_of_coeff[-1], K)

    # gcd of numerators in the coefficients
    for i in num:
        gcd_numer = i.gcd(gcd_numer)

    gcd_numer = K.new(gcd_numer.rep)

    # divide all the coefficients by the gcd
    for i, j in enumerate(list_of_coeff):
        list_of_coeff[i] = base.from_FractionField(j / gcd_numer, K)

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
    sol.append(DMFdiff(listofpoly[0]))

    for i, j in enumerate(listofpoly[1:]):
        sol.append(DMFdiff(j) + listofpoly[i])

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

    >>> from sympy.holonomic.holonomic import HolonomicFunction, DifferentialOperators
    >>> from sympy.polys.domains import ZZ, QQ
    >>> from sympy import symbols
    >>> x = symbols('x')
    >>> R, Dx = DifferentialOperators(QQ.old_poly_ring(x),'Dx')

    >>> p = HolonomicFunction(Dx - 1, x, 0, [1])  # e^x
    >>> q = HolonomicFunction(Dx**2 + 1, x, 0, [0, 1])  # sin(x)

    >>> p + q  # annihilator of e^x + sin(x)
    HolonomicFunction((-1) + (1)Dx + (-1)Dx**2 + (1)Dx**3, x), f(0) = 1 , f'(0) = 2 , f''(0) = 1

    >>> p * q  # annihilator of e^x * sin(x)
    HolonomicFunction((2) + (-2)Dx + (1)Dx**2, x), f(0) = 0 , f'(0) = 1
    """

    _op_priority = 20

    def __init__(self, annihilator, x, x0=0, y0=[]):
        """
        Takes the annihilator and variable of the function.
        x0 is the point for which initial conditions are given and
        y0 is a vector of initial values y0 = [f(x0), f'(x0), f''(x0) ...]
        To make the function unique, length of the vector `y0` must be equal to or
        greater than the order of differential equation.
        """
        # initial conditions as a list [f(x0), f'(x0), ...]
        if not isinstance(y0, list):
            self.y0 = [y0]
        else:
            self.y0 = y0

        if len(self.y0) == 0:
            self._have_init_cond = False
        else:
            self._have_init_cond = True
        # the point for initial conditions, defualt is zero.
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
                cond_str += ', f%s(%s) = %s ' % (diff_str, sstr(self.x0), sstr(i))
                diff_str += "'"

            sol = str_sol + cond_str
            return sol

    __str__ = __repr__

    def __add__(self, other):
        deg1 = self.annihilator.order
        deg2 = other.annihilator.order
        dim = max(deg1, deg2)
        R = self.annihilator.parent.base
        K = R.get_field()

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
                    p.append(K.new(expr.listofpoly[i].rep))
            r.append(p)

        r = NewMatrix(r).transpose()

        homosys = [[S(0) for q in range(dim + 1)]]
        homosys = NewMatrix(homosys).transpose()

        # solving the linear system using gauss jordan solver
        solcomp = r.gauss_jordan_solve(homosys)
        sol = solcomp[0]

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
                        p.append(S(0))
                    else:

                        p.append(K.new(expr.listofpoly[i].rep))
                r.append(p)

            r = NewMatrix(r).transpose()

            homosys = [[S(0) for q in range(dim + 1)]]
            homosys = NewMatrix(homosys).transpose()

            solcomp = r.gauss_jordan_solve(homosys)
            sol = solcomp[0]

        # taking only the coefficients needed to multiply with `self`
        # can be also be done the other way by taking R.H.S and multiplying with
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
                # initial conditions for different points
                # to be implemented
                pass

        return HolonomicFunction(sol, self.x)

    def integrate(self, *args):
        """
        Integrate the given holonomic function. Limits can be provided,
        Initial conditions can only be computed when limits are (x0, x).

        Examples
        ========

        >>> from sympy.holonomic.holonomic import HolonomicFunction, DifferentialOperators
        >>> from sympy.polys.domains import ZZ, QQ
        >>> from sympy import symbols
        >>> x = symbols('x')
        >>> R, Dx = DifferentialOperators(QQ.old_poly_ring(x),'Dx')

        >>> HolonomicFunction(Dx - 1, x, 0, 1).integrate((x, 0, x))  # e^x - 1
        HolonomicFunction((-1)Dx + (1)Dx**2, x), f(0) = 0 , f'(0) = 1

        # integrate(cos(x), (x 0, x)) = sin(x)
        >>> HolonomicFunction(Dx**2 + 1, x, 0, [1, 0]).integrate((x, 0, x))
        HolonomicFunction((1)Dx + (1)Dx**3, x), f(0) = 0 , f'(0) = 1 , f''(0) = 0
        """

        # just multiply by Dx from right
        D = self.annihilator.parent.derivative_operator
        if (not args) or (not self._have_init_cond):
            return HolonomicFunction(self.annihilator * D, self.x)

        # definite integral if limits are (x0, x)
        if len(args) == 1 and len(args[0]) == 3 and args[0][1] == self.x0 and args[0][2] == self.x:
            y0 = [S(0)]
            y0 += self.y0
            return HolonomicFunction(self.annihilator * D, self.x, self.x0, y0)

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

        if not isinstance(other, HolonomicFunction):
            if not self._have_init_cond:
                return self
            else:
                y0 = _extend_y0(self, ann_self.order)
                y1 = []
                for j in y0:
                    y1.append(j * other)
                return HolonomicFunction(ann_self, self.x, self.x0, y1)

        ann_other = other.annihilator
        list_self = []
        list_other = []
        a = ann_self.order
        b = ann_other.order

        R = ann_self.parent.base
        K = R.get_field()

        for j in ann_self.listofpoly:
            list_self.append(K.new(j.rep))

        for j in ann_other.listofpoly:
            list_other.append(K.new(j.rep))

        # will be used to reduce the degree
        self_red = [-list_self[i] / list_self[a] for i in range(a)]

        other_red = [-list_other[i] / list_other[b] for i in range(b)]

        # coeff_mull[i][j] is the coefficient of Dx^i(f).Dx^j(g)
        coeff_mul = [[S(0) for i in range(b + 1)] for j in range(a + 1)]
        coeff_mul[0][0] = S(1)

        # making the ansatz
        lin_sys = [[coeff_mul[i][j] for i in range(a) for j in range(b)]]

        homo_sys = [[S(0) for q in range(a * b)]]
        homo_sys = NewMatrix(homo_sys).transpose()

        sol = (NewMatrix(lin_sys).transpose()).gauss_jordan_solve(homo_sys)

        # until a non trivial solution is found
        while sol[0].is_zero:

            # updating the coefficents Dx^i(f).Dx^j(g) for next degree
            for i in range(a - 1, -1, -1):
                for j in range(b - 1, -1, -1):
                    coeff_mul[i][j + 1] += coeff_mul[i][j]
                    coeff_mul[i + 1][j] += coeff_mul[i][j]
                    if isinstance(coeff_mul[i][j], K.dtype):
                        coeff_mul[i][j] = DMFdiff(coeff_mul[i][j])
                    else:
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

            sol = (NewMatrix(lin_sys).transpose()).gauss_jordan_solve(homo_sys)


        sol_ann = _normalize(sol[0][0:], self.annihilator.parent, negative=False)

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

    __rmul__ = __mul__

    def __sub__(self, other):
        return self + other * -1

    def __rsub__(self, other):
        return self * -1 + other

    def __pow__(self, n):
        if n == 0:
            return S(1)
        if n == 1:
            return self
        else:
            if n % 2 == 1:
                powreduce = self**(n - 1)
                return powreduce * self
            elif n % 2 == 0:
                powreduce = self**(n / 2)
                return powreduce * powreduce

    def composition(self, expr, *args):
        """
        Returns the annihilator after composition of a holonomic function with
        an algebraic function. Initial conditions for the annihilator after composition
        can be also be provided to the function.

        Examples
        ========

        >>> from sympy.holonomic.holonomic import HolonomicFunction, DifferentialOperators
        >>> from sympy.polys.domains import ZZ, QQ
        >>> from sympy import symbols
        >>> x = symbols('x')
        >>> R, Dx = DifferentialOperators(QQ.old_poly_ring(x),'Dx')

        >>> HolonomicFunction(Dx - 1, x).composition(x**2, 0, [1])  # e^(x**2)
        HolonomicFunction((-2*x) + (1)Dx, x), f(0) = 1

        >>> HolonomicFunction(Dx**2 + 1, x).composition(x**2 - 1, 1, [1, 0])
        HolonomicFunction((4*x**3) + (-1)Dx + (x)Dx**2, x), f(1) = 1 , f'(1) = 0

        See Also
        ========

        from_hyper
        """

        R = self.annihilator.parent
        a = self.annihilator.order
        diff = expr.diff()
        listofpoly = self.annihilator.listofpoly

        for i, j in enumerate(listofpoly):
            if isinstance(j, self.annihilator.parent.base.dtype):
                listofpoly[i] = self.annihilator.parent.base.to_sympy(j)

        r = listofpoly[a].subs({self.x:expr})
        subs = [-listofpoly[i].subs({self.x:expr}) / r for i in range (a)]
        coeffs = [S(0) for i in range(a)]  # coeffs[i] == coeff of (D^i f)(a) in D^k (f(a))
        coeffs[0] = S(1)
        system = [coeffs]
        homogeneous = Matrix([[S(0) for i in range(a)]]).transpose()
        sol = S(0)

        while sol.is_zero:
            coeffs_next = [p.diff() for p in coeffs]
            for i in range(a - 1):
                coeffs_next[i + 1] += (coeffs[i] * diff)

            for i in range(a):
                coeffs_next[i] += (coeffs[-1] * subs[i] * diff)
            coeffs = coeffs_next

            # check for linear relations
            system.append(coeffs)
            sol_tuple = (Matrix(system).transpose()).gauss_jordan_solve(homogeneous)
            sol = sol_tuple[0]

        tau = sol.atoms(Dummy).pop()
        sol = sol.subs(tau, 1)
        sol = _normalize(sol[0:], R, negative=False)

        # if initial conditions are given for the resulting function
        if args:
            return HolonomicFunction(sol, self.x, args[0], args[1])
        return HolonomicFunction(sol, self.x)

    def to_sequence(self):
        """
        Finds the recurrence relation in power series expansion
        of the function.

        Examples
        ========

        >>> from sympy.holonomic.holonomic import HolonomicFunction, DifferentialOperators
        >>> from sympy.polys.domains import ZZ, QQ
        >>> from sympy import symbols
        >>> x = symbols('x')
        >>> R, Dx = DifferentialOperators(QQ.old_poly_ring(x),'Dx')

        >>> HolonomicFunction(Dx - 1, x, 0, [1]).to_sequence()
        HolonomicSequence((-1) + (n + 1)Sn, n), u(0) = 1

        See Also
        ========

        HolonomicFunction.series

        References
        ==========

        hal.inria.fr/inria-00070025/document
        """

        dict1 = {}
        n = symbols('n', integer=True)
        dom = self.annihilator.parent.base.dom
        R, _ = RecurrenceOperators(dom.old_poly_ring(n), 'Sn')

        for i, j in enumerate(self.annihilator.listofpoly):
            listofdmp = j.all_coeffs()
            degree = len(listofdmp) - 1
            for k in range(degree + 1):
                coeff = listofdmp[degree - k]
                if coeff == 0:
                    continue
                if i - k in dict1.keys():
                    dict1[i - k] += (coeff * rf(n - k + 1, i))
                else:
                    dict1[i - k] = (coeff * rf(n - k + 1, i))

        sol = []
        lower = min(dict1.keys())
        upper = max(dict1.keys())

        for j in range(lower, upper + 1):
            if j in dict1.keys():
                sol.append(dict1[j].subs(n, n - lower))
            else:
                sol.append(S(0))
        # recurrence relation
        sol = RecurrenceOperator(sol, R)

        if not self._have_init_cond:
            return HolonomicSequence(sol)
        if self.x0 != 0:
            return HolonomicSequence(sol)
        # computing the initial conditions for recurrence
        order = sol.order
        all_roots = roots(sol.listofpoly[-1].rep, filter='Z')
        all_roots = all_roots.keys()

        if all_roots:
            max_root = max(all_roots)
            if max_root >= 0:
                order += max_root + 1

        y0 = _extend_y0(self, order)
        u0 = []
        # u(n) = y^n(0)/factorial(n)
        for i, j in enumerate(y0):
            u0.append(j / factorial(i))

        return HolonomicSequence(sol, u0)

    def series(self, n=6, coefficient=False, order=True):
        """
        Finds the power series expansion of given holonomic function.

        Examples
        ========

        >>> from sympy.holonomic.holonomic import HolonomicFunction, DifferentialOperators
        >>> from sympy.polys.domains import ZZ, QQ
        >>> from sympy import symbols
        >>> x = symbols('x')
        >>> R, Dx = DifferentialOperators(QQ.old_poly_ring(x),'Dx')

        >>> HolonomicFunction(Dx - 1, x, 0, [1]).series()  # e^x
        1 + x + x**2/2 + x**3/6 + x**4/24 + x**5/120 + O(x**6)

        >>> HolonomicFunction(Dx**2 + 1, x, 0, [0, 1]).series(n=8)  # sin(x)
        x - x**3/6 + x**5/120 - x**7/5040 + O(x**8)

        See Also
        ========

        HolonomicFunction.to_sequence
        """

        recurrence = self.to_sequence()
        l = len(recurrence.u0) - 1
        k = recurrence.recurrence.order
        x = self.x
        seq_dmp = recurrence.recurrence.listofpoly
        R = recurrence.recurrence.parent.base
        K = R.get_field()
        seq = []

        if 0 in roots(seq_dmp[-1].rep, filter='Z').keys():
            singular = True
        else:
            singular = False

        for i, j in enumerate(seq_dmp):
            seq.append(K.new(j.rep))

        sub = [-seq[i] / seq[k] for i in range(k)]
        sol = [i for i in recurrence.u0]

        if l + 1 >= n:
            pass
        else:
            # use the initial conditions to find the next term
            for i in range(l + 1 - k, n - k):
                coeff = S(0)
                for j in range(k):
                    if i + j >= 0:
                        coeff += DMFsubs(sub[j], i) * sol[i + j]
                sol.append(coeff)

        if coefficient:
            return sol

        ser = S(0)
        for i, j in enumerate(sol):
            ser += x**i * j
        if order:
            return ser + Order(x**n, x)
        else:
            return ser

    def _indicial(self):
        """Computes the roots of Indicial equation.
        """

        list_coeff = self.annihilator.listofpoly
        R = self.annihilator.parent.base
        x = self.x
        s = R.zero
        y = R.one

        def _pole_degree(poly):
            root_all = roots(poly.rep, filter='Z')
            if 0 in root_all.keys():
                return root_all[0]
            else:
                return 0

        degree = [j.degree() for j in list_coeff]
        degree = max(degree)
        inf = 10 * (max(1, degree) + max(1, self.annihilator.order))

        deg = lambda q: inf if q.is_zero else _pole_degree(q)
        b = deg(list_coeff[0])
        print (b)

        for j in range(1, len(list_coeff)):
            b = min(b, deg(list_coeff[j]) - j)
            print(b)

        for i, j in enumerate(list_coeff):
            listofdmp = j.all_coeffs()
            degree = len(listofdmp) - 1
            if - i - b <= 0:
                s = s + listofdmp[degree - i - b] * y
            y *= x - i
        return roots(s.rep, filter='R').keys()

    def evalf(self, points, method='RK4'):
        """
        Finds numerical value of a holonomic function using numerical methods.
        (RK4 by default). A set of points (real or complex) must be provided
        which will be the path for the numerical integration.

        The path should be given as a list [x1, x2, ... xn]. The numerical values
        will be computed at each point in this order x1 --> x2 --> x3 ... --> xn.

        Returns values of the function at x1, x2, ... xn in a list.

        Examples
        =======

        >>> from sympy.holonomic.holonomic import HolonomicFunction, DifferentialOperators
        >>> from sympy.polys.domains import ZZ, QQ
        >>> from sympy import symbols
        >>> x = symbols('x')
        >>> R, Dx = DifferentialOperators(QQ.old_poly_ring(x),'Dx')
        >>> # a straight line on the real axis from (0 to 1)
        >>> r = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]

        # using Runge-Kutta 4th order on e^x from 0.1 to 1.
        # exact solution at 1 is 2.71828182845905
        >>> HolonomicFunction(Dx - 1, x, 0, [1]).evalf(r)
        [1.10517083333333, 1.22140257085069, 1.34985849706254, 1.49182424008069, 
        1.64872063859684, 1.82211796209193, 2.01375162659678, 2.22553956329232, 
        2.45960141378007, 2.71827974413517]

        # using Euler's method for the same
        >>> HolonomicFunction(Dx - 1, x, 0, [1]).evalf(r, method='Euler')
        [1.1, 1.21, 1.331, 1.4641, 1.61051, 1.771561, 1.9487171, 2.14358881,
         2.357947691, 2.5937424601]

        One can also observe that the value obtained using Runge-Kutta 4th order
        is much more accurate than Euler's method.
        """

        from sympy.holonomic.numerical import _evalf
        return _evalf(self, points, method=method)


def from_hyper(func, x0=0, evalf=False):
    """
    Converts Hypergeometric Function to Holonomic.
    func is the Hypergeometric Function and x0 be the point at
    which initial conditions are required.

    Examples
    =======

    >>> from sympy.holonomic.holonomic import from_hyper, DifferentialOperators
    >>> from sympy import symbols, hyper, S
    >>> x = symbols('x')
    >>> from_hyper(hyper([], [S(3)/2], x**2/4))
    HolonomicFunction((-x) + (2)Dx + (x)Dx**2, x), f(1) = sinh(1) , f'(1) = -sinh(1) + cosh(1)
    """

    a = func.ap
    b = func.bq
    z = func.args[2]
    x = z.atoms(Symbol).pop()
    R, Dx = DifferentialOperators(QQ.old_poly_ring(x), 'Dx')

    # generalized hypergeometric differential equation
    r1 = 1
    for i in range(len(a)):
        r1 = r1 * (x * Dx + a[i])
    r2 = Dx
    for i in range(len(b)):
        r2 = r2 * (x * Dx + b[i] - 1)
    sol = r1 - r2

    simp = hyperexpand(func)

    if isinstance(simp, Infinity) or isinstance(simp, NegativeInfinity):
        return HolonomicFunction(sol, x).composition(z)

    def _find_conditions(simp, x, x0, order, evalf=False):
        y0 = []
        for i in range(order):
            if evalf:
                val = simp.subs(x, x0).evalf()
            else:
                val = simp.subs(x, x0)
            # return None if it is Infinite or NaN
            if (val.is_finite is not None and not val.is_finite) or isinstance(val, NaN):
                return None
            y0.append(val)
            simp = simp.diff()
        return y0

    # if the function is known symbolically
    if not isinstance(simp, hyper):
        y0 = _find_conditions(simp, x, x0, sol.order)
        while not y0:
            # if values don't exist at 0, then try to find initial
            # conditions at 1. If it doesn't exist at 1 too then
            # try 2 and so on.
            x0 += 1
            y0 = _find_conditions(simp, x, x0, sol.order)

        return HolonomicFunction(sol, x).composition(z, x0, y0)

    if isinstance(simp, hyper):
        x0 = 1
        # use evalf if the function can't be simpified
        y0 = _find_conditions(simp, x, x0, sol.order, evalf)
        while not y0:
            x0 += 1
            y0 = _find_conditions(simp, x, x0, sol.order, evalf)
        return HolonomicFunction(sol, x).composition(z, x0, y0)

    return HolonomicFunction(sol, x).composition(z)


def from_meijerg(func, x0=0, evalf=False):
    """
    Converts a Meijer G-function to Holonomic.
    func is the Hypergeometric Function and x0 be the point at
    which initial conditions are required.

    Examples
    =======

    >>> from sympy.holonomic.holonomic import from_meijerg, DifferentialOperators
    >>> from sympy import symbols, meijerg, S
    >>> x = symbols('x')
    >>> from_meijerg(meijerg(([], []), ([S(1)/2], [0]), x**2/4))
    HolonomicFunction((1) + (1)Dx**2, x), f(0) = 0 , f'(0) = 1/sqrt(pi)
    """

    a = func.ap
    b = func.bq
    n = len(func.an)
    m = len(func.bm)
    p = len(a)
    z = func.args[2]
    x = z.atoms(Symbol).pop()
    R, Dx = DifferentialOperators(QQ.old_poly_ring(x), 'Dx')

    # compute the differential equation satisfied by the
    # Meijer G-function.
    mnp = (-1)**(m + n - p)
    r1 = x * mnp

    for i in range(len(a)):
        r1 *= x * Dx + 1 - a[i]

    r2 = 1

    for i in range(len(b)):
        r2 *= x * Dx - b[i]

    sol = r1 - r2

    simp = hyperexpand(func)

    if isinstance(simp, Infinity) or isinstance(simp, NegativeInfinity):
        return HolonomicFunction(sol, x).composition(z)

    def _find_conditions(simp, x, x0, order, evalf=False):
        y0 = []
        for i in range(order):
            if evalf:
                val = simp.subs(x, x0).evalf()
            else:
                val = simp.subs(x, x0)
            if (val.is_finite is not None and not val.is_finite) or isinstance(val, NaN):
                return None
            y0.append(val)
            simp = simp.diff()
        return y0

    # computing initial conditions
    if not isinstance(simp, meijerg):
        y0 = _find_conditions(simp, x, x0, sol.order)
        while not y0:
            x0 += 1
            y0 = _find_conditions(simp, x, x0, sol.order)

        return HolonomicFunction(sol, x).composition(z, x0, y0)

    if isinstance(simp, meijerg):
        x0 = 1
        y0 = _find_conditions(simp, x, x0, sol.order, evalf)
        while not y0:
            x0 += 1
            y0 = _find_conditions(simp, x, x0, sol.order, evalf)

        return HolonomicFunction(sol, x).composition(z, x0, y0)

    return HolonomicFunction(sol, x).composition(z)


def _extend_y0(Holonomic, n):
    """
    Tries to find more initial conditions by substituting the initial
    value point in the differential equation.
    """

    annihilator = Holonomic.annihilator
    a = annihilator.order
    x = Holonomic.x
    listofpoly = []
    y0 = Holonomic.y0
    R = annihilator.parent.base
    K = R.get_field()

    for i, j in enumerate(annihilator.listofpoly):
            if isinstance(j, annihilator.parent.base.dtype):
                listofpoly.append(K.new(j.rep))

    if len(y0) < a or n <= len(y0):
        return y0
    else:
        list_red = [-listofpoly[i] / listofpoly[a]
                    for i in range(a)]
        y1 = [i for i  in y0]
        for i in range(n - a):
            sol = 0
            for a, b in zip(y1, list_red):
                r = DMFsubs(b, Holonomic.x0)
                if not r.is_finite:
                    return y0
                sol += a * r
            y1.append(sol)
            list_red = _derivate_diff_eq(list_red)

        return y0 + y1[len(y0):]


def DMFdiff(frac):
    # differentiate a p/q DMF object
    if not isinstance(frac, DMF):
        return frac.diff()

    K = frac.ring
    p = K.numer(frac)
    q = K.denom(frac)
    sol_num = - p * q.diff() + q * p.diff()
    sol_denom = q**2
    return K((sol_num.rep, sol_denom.rep))


def DMFsubs(frac, x0):
    # substitute the point x0 in DMF object of the form p/q
    if not isinstance(frac, DMF):
        return frac

    p = frac.num
    q = frac.den
    sol_p = S(0)
    sol_q = S(0)

    for i, j in enumerate(reversed(p)):
        if isinstance(j, PythonRational):
            j = sympify(j)
        sol_p += j * x0**i

    for i, j in enumerate(reversed(q)):
        if isinstance(j, PythonRational):
            j = sympify(j)
        sol_q += j * x0**i

    return sol_p / sol_q


def from_sympy(func):
    """
    Uses `meijerint._rewrite1` to convert to `meijerg` function and then eventually
    to Holonomic Functions. Only works when `meijerint._rewrite1` returns a `meijerg`
    representation of the function provided.

    Examples
    ========

    >>> from sympy.holonomic.holonomic import from_sympy
    >>> from sympy import sin, exp, symbols
    >>> x = symbols('x')
    >>> from_sympy(sin(x))
    HolonomicFunction((1) + (1)Dx**2, x), f(0) = 0 , f'(0) = 1

    >>> from_sympy(exp(x))
    HolonomicFunction((-1) + (1)Dx, x), f(0) = 1

    See Also
    ========

    meijerint._rewrite1
    """

    x = func.atoms(Symbol).pop()
    args = meijerint._rewrite1(func, x)

    if args:
        fac, po, g, _ = args
    else:
        return None

    # lists for sum of meijerg functions
    fac_list = [fac * i[0] for i in g]
    t = po.as_base_exp()
    s = t[1] if t[0] is x else S(0)
    po_list = [s + i[1] for i in g]
    G_list = [i[2] for i in g]

    # finds meijerg representation of x**s * meijerg(a1 ... ap, b1 ... bq, z)
    def _shift(func, s):
        z = func.args[-1]
        d = z.collect(x, evaluate=False)
        b = list(d)[0]
        a = d[b]

        if isinstance(a, exp_polar):
            a = exp(a.as_base_exp()[1])
            z = a * b

        t = b.as_base_exp()
        b = t[1] if t[0] is x else S(0)
        r = s / b
        an = (i + r for i in func.args[0][0])
        ap = (i + r for i in func.args[0][1])
        bm = (i + r for i in func.args[1][0])
        bq = (i + r for i in func.args[1][1])

        return a**-r, meijerg((an, ap), (bm, bq), z)

    coeff, m = _shift(G_list[0], po_list[0])
    sol = fac_list[0] * coeff * from_meijerg(m)

    # add all the meijerg functions after converting to holonomic
    for i in range(1, len(G_list)):
        coeff, m = _shift(G_list[i], po_list[i])
        sol += fac_list[i] * coeff * from_meijerg(m)

    return sol
