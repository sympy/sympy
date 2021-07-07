r"""
This module contains :py:meth:`~sympy.solvers.ode.riccati.solve_riccati`,
a function which gives all rational particular solutions to first order
Riccati ODEs. A general first order Riccati ODE is given by -

.. math:: y' = b_0(x) + b_1(x)w + b_2(x)w^2

where `b_0, b_1` and `b_2` can be arbitrary rational functions of `x`
with `b_2 \ne 0`.

Background
==========

A Riccati equation can be transformed to its normal form

.. math:: y' + y^2 = a(x)

using the transformation

.. math:: y = -b_2(x) - \frac{b'_2(x)}{2 b_2(x)} - \frac{b_1(x)}{2}

where `a(x)` is given by

.. math:: a(x) = \frac{1}{4}\left(\frac{b_2'}{b_2} + b_1\right)^2 - \frac{1}{2}\left(\frac{b_2'}{b_2} + b_1\right)' - b_0 b_2

Thus, we can develop an algorithm to solve for the Riccati equation
in its normal form, which would in turn give us the solution for
the original Riccati equation.

Algorithm
=========

The algorithm implemented here is presented in the Ph.D thesis
"Rational and Algebraic Solutions of First-Order Algebraic ODEs"
by N. Thieu Vo. The entire thesis can be found here -
https://www3.risc.jku.at/publications/download/risc_5387/PhDThesisThieu.pdf

We have only implemented the Rational Riccati solver (Algorithm 11,
Pg 78-82 in Thesis). Before we proceed towards the implementation
of the algorithm, a few definitions to understand are -

1. Valuation of a Rational Function at `\infty`:
    The valuation of a rational function `p(x)` at `\infty` is equal
    to the difference between the degree of the denominator and the
    numerator of `p(x)`.

    NOTE: A general definition of valuation of a rational function
    at any value of `x` can be found in Pg 63 of the thesis, but
    is not of any interest for this algorithm.

2. Zeros and Poles of a Rational Function:
    Let `a(x) = \frac{S(x)}{T(x)}, T \ne 0` be a rational function
    of `x`. Then -

    a. The Zeros of `a(x)` are the roots of `S(x)`.
    b. The Poles of `a(x)` are the roots of `T(x)`. However, `\infty`
    can also be a pole of a(x). We say that `a(x)` has a pole at
    `\infty` if `a(\frac{1}{x})` has a pole at 0.

Every pole is associated with an order that is equal to the multiplicity
of its appearence as a root of `T(x)`. A pole is called a simple pole if
it has an order 1. Similarly, a pole is called a multiple pole if it has
an order `\ge` 2.

Necessary Conditions
====================

For a Riccati equation in its normal form,

.. math:: y' + y^2 = a(x)

we can define

a. A pole is called a movable pole if it is a pole of `y(x)` and is not
a pole of `a(x)`.
b. Similarly, a pole is called a non-movable pole if it is a pole of both
`y(x)` and `a(x)`.

Then, the algorithm states that a rational solution exists only if -

a. Every pole of `a(x)` must be either a simple pole or a multiple pole
of even order.
b. The valuation of `a(x)` at `\infty` must be even or be `\ge` 2.

This algorithm finds all possible rational solutions for the Riccati ODE.
If no rational solutions are found, it means that no rational solutions
exist.

Solution
========

With these definitions, we can state a general form for the solution of
the equation. `y(x)` must have the form -

.. math:: y(x) = \sum_{i=1}^{n} \sum_{j=1}^{r_i} \frac{c_{ij}}{(x - x_i)^j} + \sum_{i=1}^{m} \frac{1}{x - \chi_i} + \sum_{i=0}^{N} d_i x^i

where `x_1, x_2, ..., x_n` are non-movable poles of `a(x)`,
`\chi_1, \chi_2, ..., \chi_m` are movable poles of `a(x)`, and the values
of `N, n, r_1, r_2, ..., r_n` can be determined from `a(x)`. The
coefficient vectors `(d_0, d_1, ..., d_N)` and `(c_{i1}, c_{i2}, ..., c_{i r_i})`
can be determined from `a(x)`. We will have 2 choices each of these vectors
and part of the procedure is figuring out which of the 2 should be used
to get the solution correctly.

Implementation
==============

The code is written to match the steps given in the thesis (Pg 82)

Step 0 : Match the equation -
Find `b_0, b_1` and `b_2`. If `b_2 = 0` or no such functions exist, raise
an error

Step 1 : Transform the equation to its normal form as explained in the
theory section.

Step 2 : Initialize 2 empty set of solutions, ``sol`` and ``presol``.

Step 3 : If `a(x) = 0`, append `\frac{1}/{(x - C1)}` to ``presol``.

Step 4 : If `a(x)` is a rational non-zero number, append `\pm \sqrt{a}`
to ``presol``.

Step 5 : Find the poles and their multiplicities of `a(x)` using
``find_poles``. Let the number of poles be `n`. Also find the
valuation of `a(x)` at `\infty` using ``val_at_inf``.

NOTE: Although the algorithm considers `\infty` as a pole, it is
not mentioned if it a part of the set of finite poles. `\infty`
is NOT a part of the set of finite poles. If a pole exists at
`\infty`, we use its multiplicty to find the laurent series of
`a(x)` about `\infty`.

Step 6 : Find `n` c-vectors (one for each pole) and 1 d-vector using
``construct_c`` and ``construct_d``. Now, determine all the ``2**(n + 1)``
combinations of choosing between 2 choices for each of the `n` c-vectors
and 1 d-vector.

NOTE: The equation for `d_{-1}` in Case 4 (Pg 80) has a printinig
mistake. The term `- d_N` must be replaced with `-N d_N`. The same
has been explained in the code as well.

For each of these above combinations, do

Step 8 : Compute `m` in ``compute_m_ybar``

Step 9 : In ``compute_m_ybar``, compute ybar as well.

Step 10 : If `m` is a non-negative integer -

Step 11: Find a polynomial solution of degree `m` for the auxiliary equation.

There are 2 cases possible -

    a. `m` is a non-negative integer: We can solve for the coefficients
    in `p(x)` using Undetermined Coefficients.

    b. `m` is an expression: In this case, we may not be able to determine
    if it is possible for `m` to be a non-negative integer. However, the
    solution to the auxiliary equation always seems to be x**m + C1 (where
    C1 is an arbitrary constant). This has been confirmed by testing many
    examples. NOTE that this is NOT mentioned in the thesis and is an
    assumption made by me. If any failing test cases are found, this case
    must be removed.

Step 12 : For each `p(x)` that exists, append `ybar + \frac{p'(x)}{p(x)}`
to ``presol``.

Step 13 : For each solution in ``presol``, apply an inverse transformation
and append them to ``sol``, so that the solutions of the original equation
are found using the solutions of the equation in its normal form.
"""

from itertools import product
from sympy.core import S
from sympy.core.add import Add
from sympy.core.numbers import oo
from sympy.core.function import count_ops
from sympy.core.relational import Eq
from sympy.core.symbol import symbols, Symbol, Dummy
from sympy.functions import sqrt, exp
from sympy.functions.elementary.complexes import sign
from sympy.integrals.integrals import Integral
from sympy.polys.domains import ZZ
from sympy.polys.polytools import Poly
from sympy.polys.polyroots import roots
from sympy.solvers.solveset import linsolve
from sympy.tensor.indexed import IndexedBase, Indexed

def riccati_normal(w, x, b1, b2):
    # y(x) = -b2(x)*w(x) - b2'(x)/(2*b2(x)) - b1(x)/2
    return -b2*w - b2.diff(x)/(2*b2) - b1/2


def riccati_inverse_normal(y, x, b1, b2, bp=None):
    # bp is the expression which is independent of the solution
    # and hence, it need not be computed again
    if bp is None:
        bp = -b2.diff(x)/(2*b2**2) - b1/(2*b2)
    # w(x) = -y(x)/b2(x) - b2'(x)/(2*b2(x)^2) - b1(x)/(2*b2(x))
    return -y/b2 + bp


def riccati_reduced(eq, f, x):
    match, funcs = match_riccati(eq, f, x)
    # If equation is not a Riccati ODE, exit
    if not match:
        return False
    # Using the rational functions, find the expression for a(x)
    b0, b1, b2 = funcs
    a = -b0*b2 + b1**2/4 - b1.diff(x)/2 + 3*b2.diff(x)**2/(4*b2**2) + b1*b2.diff(x)/(2*b2) - \
        b2.diff(x, 2)/(2*b2)
    # Normal form of Riccati ODE is f'(x) + f(x)^2 = a(x)
    return f(x).diff(x) + f(x)**2 - a

def linsolve_dict(eq, syms):
    # Convert tuple type return value of linsolve
    # to a dictionary for ease of use
    sol = linsolve(eq, syms)
    if not sol:
        return {}
    return {k:v for k, v in zip(syms, list(sol)[0])}


def match_riccati(eq, f, x):
    # Group terms based on f(x)
    eq = eq.expand().collect(f(x))
    cf = eq.coeff(f(x).diff(x))

    # There must be an f(x).diff(x) term.
    # eq must be an Add object since we are using the expanded
    # equation and it must have atleast 2 terms (b2 != 0)
    if cf != 0 and isinstance(eq, Add):

        # Divide all coefficients by the coefficient of f(x).diff(x)
        # and add the terms again to get the same equation
        eq = Add(*map(lambda x: (x/cf).cancel(), eq.args)).collect(f(x))

        # Match the equation with the pattern
        b1 = -eq.coeff(f(x))
        b2 = -eq.coeff(f(x)**2)
        b0 = (f(x).diff(x) - b1*f(x) - b2*f(x)**2 - eq).expand()

        # If b_0(x) contains f(x), it is not a Riccati ODE
        if len(b0.atoms(f)) or not all([b2 != 0, b0.is_rational_function(x), \
            b1.is_rational_function(x), b2.is_rational_function(x)]):
            return False, []
        return True, [b0, b1, b2]
    return False, []


def find_poles(den, x):
    # All roots of denominator are poles of a(x)
    p = roots(den, x)
    return p


def val_at_inf(num, den, x):
    # Valuation of a rational function at oo = deg(denom) - deg(numer)
    return den.degree(x) - num.degree(x)


def check_necessary_conds(val_inf, muls):
    """
    The necessary conditions for a rational solution
    to exist are as follows -

    i) Every pole of a(x) must be either a simple pole
    or a multiple pole of even order.

    ii) The valuation of a(x) at infinity must be even
    or be greater than or equal to 2.

    Here, a simple pole is a pole with multiplicity 1
    and a multiple pole is a pole with multiplicity
    greater than 1.
    """
    return (val_inf >= 2 or (val_inf <= 0 and val_inf%2 == 0)) and \
        all([mul == 1 or (mul%2 == 0 and mul >= 2) for mul in muls])


def inverse_transform_poly(num, den, x):
    """
    A function to make the substitution
    x -> 1/x in Poly objects.
    """
    # Declare for reuse
    one = Poly(1, x)
    xpoly = Poly(x, x)

    # Check if degree of numerator is same as denominator
    pwr = val_at_inf(num, den, x)
    if pwr >= 0:
        # Denominator has greater degree. Substituting x with
        # 1/x would make the extra power go to the numerator
        if num.expr != 0:
            num = num.transform(one, xpoly) * x**pwr
            den = den.transform(one, xpoly)
    else:
        # Numerator has greater degree. Substituting x with
        # 1/x would make the extra power go to the denominator
        num = num.transform(one, xpoly)
        den = den.transform(one, xpoly) * x**(-pwr)
    return num.cancel(den, include=True)


def limit_at_inf(num, den, x):
    # pwr = degree(num) - degree(den)
    pwr = -val_at_inf(num, den, x)
    # Numerator has a greater degree than denominator
    # Limit at infinity would depend on the sign of the
    # leading coefficients of numerator and denominator
    if pwr > 0:
        return oo*sign(num.LC()/den.LC())
    # Degree of numerator is equal to that of denominator
    # Limit at infinity is just the ratio of leading coeffs
    elif pwr == 0:
        return num.LC()/den.LC()
    # Degree of numerator is less than that of denominator
    # Limit at infinity is just 0
    else:
        return 0


def construct_c_case_1(num, den, x, pole):
    # Find the coefficient of 1/(x - pole)**2 in the
    # Laurent series expansion of a(x) about pole.
    num1, den1 = (num*Poly((x - pole)**2, x, extension=True)).cancel(den, include=True)
    r = (num1.subs(x, pole))/(den1.subs(x, pole))

    # If multiplicity is 2, the coefficient to be added
    # in the c-vector is c = (1 +- sqrt(1 + 4*r))/2
    return [[(1 + sqrt(1 + 4*r))/2], [(1 - sqrt(1 + 4*r))/2]]


def construct_c_case_2(num, den, x, pole, mul):
    # Generate the coefficients using the recurrence
    # relation mentioned in (5.14) in the thesis (Pg 80)

    # r_i = mul/2
    ri = mul//2

    # Find the Laurent series coefficients about the pole
    ser = rational_laurent_series(num, den, x, pole, mul, 6)

    # Start with an empty memo to store the coefficients
    # This is for the plus case
    cplus = [0 for i in range(ri)]

    # Base Case
    cplus[ri-1] = sqrt(ser[2*ri])

    # Iterate backwards to find all coefficients
    s = ri - 1
    sm = 0
    for s in range(ri-1, 0, -1):
        sm = 0
        for j in range(s+1, ri):
            sm += cplus[j - 1]*cplus[ri + s - j -1]
        if s!= 1:
            cplus[s-1] = (ser[ri + s] - sm)/(2*cplus[ri - 1])

    # Memo for the minus case
    cminus = [-x for x in cplus]

    # Find the 0th coefficient in the recurrence
    cplus[0] = (ser[ri + s] - sm - ri*cplus[ri - 1])/(2*cplus[ri - 1])
    cminus[0] = (ser[ri + s] - sm  - ri*cminus[ri - 1])/(2*cminus[ri - 1])

    # Add both the plus and minus cases' coefficients
    return [cplus, cminus]


def construct_c_case_3():
    # If multiplicity is 1, the coefficient to be added
    # in the c-vector is 1 (no choice)
    return [[1], [1]]


def construct_c(num, den, x, poles, muls):
    c = []
    for pole, mul in zip(poles, muls):
        c.append([])

        # Case 3
        if mul == 1:
            # Add the coefficients from Case 3
            c[-1].extend(construct_c_case_3())

        # Case 1
        elif mul == 2:
            # Add the coefficients from Case 1
            c[-1].extend(construct_c_case_1(num, den, x, pole))

        # Case 2
        else:
            # Add the coefficients from Case 2
            c[-1].extend(construct_c_case_2(num, den, x, pole, mul))

    return c


def construct_d_case_4(ser, N):
    # Initialize an empty vector
    dplus = [0 for i in range(N+2)]
    # d_N = sqrt(a_{2*N})
    dplus[N] = sqrt(ser[2*N])

    # Use the recurrence relations to find
    # the value of d_s
    for s in range(N-1, -2, -1):
        sm = 0
        for j in range(s+1, N):
            sm += dplus[j]*dplus[N+s-j]
        if s != -1:
            dplus[s] = (ser[N+s] - sm)/(2*dplus[N])

    # Coefficients for the case of d_N = -sqrt(a_{2*N})
    dminus = [-x for x in dplus]

    # The third equation in Eq 5.15 of the thesis is WRONG!
    # d_N must be replaced with N*d_N in that equation.
    dplus[-1] = (ser[N+s] - N*dplus[N] - sm)/(2*dplus[N])
    dminus[-1] = (ser[N+s] - N*dminus[N] - sm)/(2*dminus[N])

    return [dplus, dminus]


def construct_d_case_5(ser):
    # List to store coefficients for plus case
    dplus = [0, 0]

    # d_0  = sqrt(a_0)
    dplus[0] = sqrt(ser[0])

    # d_(-1) = a_(-1)/(2*d_0)
    dplus[-1] = ser[-1]/(2*dplus[0])

    # Coefficients for the minus case are just the negative
    # of the coefficients for the positive case.
    dminus = [-x for x in dplus]

    return [dplus, dminus]


def construct_d_case_6(num, den, x):
    # s_oo = lim x->0 1/x**2 * a(1/x) which is equivalent to
    # s_oo = lim x->oo x**2 * a(x)
    s_inf = limit_at_inf(Poly(x**2, x)*num, den, x)

    # d_(-1) = (1 +- sqrt(1 + 4*s_oo))/2
    return [[(1 + sqrt(1 + 4*s_inf))/2], [(1 - sqrt(1 + 4*s_inf))/2]]


def construct_d(num, den, x, val_inf):
    N = -val_inf//2
    # Multiplicity of oo as a pole
    mul = -val_inf if val_inf < 0 else 0
    ser = rational_laurent_series(num, den, x, oo, mul, 1)

    # Case 4
    if val_inf < 0:
        d = construct_d_case_4(ser, N)

    # Case 5
    elif val_inf == 0:
        d = construct_d_case_5(ser)

    # Case 6
    else:
        d = construct_d_case_6(num, den, x)

    return d


def rational_laurent_series(num, den, x, x0, mul, n):
    # Basic Idea for this function is explained here:
    # https://github.com/sympy/sympy/pull/21459#issuecomment-850500922

    # Variable for symbolic index
    m = symbols('m')
    one = Poly(1, x, extension=True)

    if x0 == oo:
        # Series at x = oo is equal to first transforming
        # the function from x -> 1/x and finding the
        # series at x = 0
        num, den = inverse_transform_poly(num, den, x)
        x0 = S(0)

    if x0:
        # For an expansion about a non-zero point, a
        # transformation from x -> x + x0 must be made
        num = num.transform(Poly(x + x0, x, extension=True), one)
        den = den.transform(Poly(x + x0, x, extension=True), one)
        num, den = num.cancel(den, include=True)

    # Remove the pole from the denominator if the series
    # expansion is about one of the poles
    num, den = (num*x**mul).cancel(den, include=True)

    # IndexedBase object to represent recurrence relation
    # for the Laurent series coefficients
    c = IndexedBase('c')
    coeff_max = max(den.degree(x), num.degree(x) + 1)
    # All coefficient values that must be found
    indices = [c[i] for i in range(-mul, max(n, coeff_max))]
    out, sums = Poly(0, x, domain=ZZ[indices + [S(x0)]]), 0

    # For each term in the denominator, add a term in
    # the recurrence relation
    for pw, cf in den.as_dict().items():
        pw = pw[0]
        for i in range(pw, coeff_max):
            out += Poly(cf*c[i - pw]*x**i, x, domain=ZZ[indices + [S(x0)]])
        sums += cf*c[m - pw]

    # If the recurrence relation is a valid expression,
    # find the base case coefficients and use these
    # to find values of other coefficients
    if out:
        coeffs = linsolve_dict((num - out).all_coeffs(), list(out.atoms(Indexed)))
        for i in range(len(coeffs), n + mul):
            ai = linsolve_dict([sums.subs(m, i).subs(coeffs)], [c[i]])
            if len(ai) > 0:
                coeffs.update({c[i]: list(ai.values())[0]})
        coeffs = list(coeffs.items())
        coeffs = {mul - k.indices[0]: v for k, v in coeffs}
    return coeffs

def compute_m_ybar(x, poles, choice, c, d, N):
    ybar = 0
    m = Poly(d[choice[0]][-1], x, extension=True)

    # Calculate the first (nested) summation for ybar
    # as given in Step 9 of the Thesis (Pg 82)
    for i in range(len(poles)):
        for j in range(len(c[i][choice[i + 1]])):
            ybar += c[i][choice[i + 1]][j]/(x - poles[i])**(j+1)
        m -= Poly(c[i][choice[i + 1]][0], x, extension=True)

    # Calculate the second summation for ybar
    for i in range(N+1):
        ybar += d[choice[0]][i]*x**i
    return (m.expr, ybar)


# def rational_laurent_series_rootsum(num, den, x):
#     pass
#     y = symbols('y')
#     coeff = lambda r, n: (num*diff((1/den).subs(x, x + r) * x, x, n) / factorial(n)).cancel(extension=True).subs(x, y)


def solve_aux_eq(numa, dena, numy, deny, x, m):
    # Assume that the solution is of the type
    # p(x) = C0 + C1*x + ... + Cm*x**m
    psyms = symbols(f'C0:{m}', cls=Dummy)
    domain = ZZ[psyms]
    psyms += (1, )
    psol = 0
    for i in range(m + 1):
        psol += Poly(psyms[i]*x**i, x, domain=domain)

    # Eq (5.16) in Thesis - Pg 81
    auxeq = (dena*(numy.diff(x)*deny - numy*deny.diff(x) + numy**2) - numa*deny**2)*psol
    if m >= 1:
        px = psol.diff(x)
        auxeq += 2*px*numy*deny*dena
    if m >= 2:
        auxeq += px.diff(x)*deny**2*dena
    if m != 0:
        # m is a non-zero integer. Find the constant terms using undetermined coefficients
        return psol, linsolve_dict(auxeq.all_coeffs(), psyms), True
    else:
        # m == 0 . Check if 1 (x**0) is a solution to the auxiliary equation
        return S(1), auxeq, auxeq == 0


def remove_redundant_sols(sol1, sol2, x):
    # If y1 and y2 are redundant solutions, there is
    # some value of the arbitrary constant for which
    # they will be equal

    syms1 = sol1.atoms(Symbol, Dummy)
    syms2 = sol2.atoms(Symbol, Dummy)
    num1, den1 = [Poly(e, x, extension=True) for e in sol1.together().as_numer_denom()]
    num2, den2 = [Poly(e, x, extension=True) for e in sol2.together().as_numer_denom()]
    # Cross multiply
    e = num1*den2 - den1*num2
    # Check if there are any constants
    syms = list(e.atoms(Symbol, Dummy))
    if len(syms):
        # Find values of constants for which solutions are equal
        redn = linsolve(e.all_coeffs(), syms)
        if len(redn):
            # Return the general solution over a particular solution
            if len(syms1) > len(syms2):
                return sol2
            # If both have constants, return the lesser complex solution
            elif len(syms1) == len(syms2):
                return sol1 if count_ops(syms1) >= count_ops(syms2) else sol2
            else:
                return sol1


def get_gen_sol_from_part_sol(part_sols, a, x):
    # There are 3 cases to find the general solution
    # from the particular solutions for a Riccati ODE
    # depending on the number of particular solution(s)
    # we have - 1, 2 or 3.
    # For more information, see Section 6 of
    # "Methods of Solution of the Riccati Differential Equation"
    # by D. R. Haaheim and F. M. Stein

    # If no particular solutions are found, a general
    # solution cannot be found
    if len(part_sols) == 0:
        return []

    # In case of a single particular solution, the general
    # solution can be found by solving a Bernoulli ODE
    elif len(part_sols) == 1:
        y1 = part_sols[0]
        z = exp(Integral(2*y1, x)) * Integral(exp(-Integral(2*y1, x)) * a, x)
        z = z.doit()
        if a == 0 or z == 0:
            return y1
        return y1 + 1/z

    # In case of 2 particular solutions, the general solution
    # can be found by solving a separable equation
    elif len(part_sols) == 2:
        y1, y2 = part_sols
        u = exp(Integral(y2 - y1, x)).doit()
        if u == 1:
            return y2
        return (y2*u - y1)/(u - 1)

    # In case of 3 particular solutions, a closed form
    # of the general solution can be obtained directly
    else:
        y1, y2, y3 = part_sols[:3]
        C1 = Dummy('C1')
        return (C1 + 1)*y2*(y1 - y3)/(C1*y1 + y2 - (C1 + 1)*y3)


def solve_riccati(fx, x, b0, b1, b2, gensol=False):
    # Step 1 : Convert to Normal Form
    a = -b0*b2 + b1**2/4 - b1.diff(x)/2 + 3*b2.diff(x)**2/(4*b2**2) + b1*b2.diff(x)/(2*b2) - \
        b2.diff(x, 2)/(2*b2)
    a_t = a.together()
    num, den = [Poly(e, x, extension=True) for e in a_t.as_numer_denom()]
    num, den = num.cancel(den, include=True)

    # Step 2
    presol = []

    # Step 3 : a(x) is 0
    if num == 0:
        presol.append(1/(x + Dummy('C1')))

    # Step 4 : a(x) is a non-zero constant
    elif x not in num.free_symbols.union(den.free_symbols):
        presol.extend([sqrt(a), -sqrt(a)])

    # Step 5 : Find poles and valuation at infinity
    poles = find_poles(den, x)
    poles, muls = list(poles.keys()), list(poles.values())
    val_inf = val_at_inf(num, den, x)

    if len(poles):
        # Check necessary conditions (outlined in the module docstring)
        if not check_necessary_conds(val_inf, muls):
            raise ValueError("Rational Solution doesn't exist")

        # Step 6
        # Construct c-vectors for each singular point
        c = construct_c(num, den, x, poles, muls)

        # Construct d vectors for each singular point
        d = construct_d(num, den, x, val_inf)

        # Step 7 : Iterate over all possible combinations and return solutions
        # For each possible combination, generate an array of 0's and 1's
        # where 0 means pick 1st choice and 1 means pick the second choice.
        choices = product(range(2), repeat=len(poles) + 1)
        m_ybars = set()
        for choice in choices:
            # Step 8 and 9 : Compute m and ybar
            m_ybars.add(compute_m_ybar(x, poles, choice, c, d, -val_inf//2))

        for m, ybar in sorted(m_ybars, key=lambda x: x[0]):
            numy, deny = [Poly(e, x, extension=True) for e in ybar.together().as_numer_denom()]
            # Step 10 : Check if a valid solution exists. If yes, also check
            # if m is a non-negative integer
            if m.is_nonnegative == True and m.is_integer == True:

                # Step 11 : Find polynomial solutions of degree m for the auxiliary equation
                psol, coeffs, exists = solve_aux_eq(num, den, numy, deny, x, m)

                # Step 12 : If valid polynomial solution exists, append solution.
                if exists:
                    # m == 0 case
                    if psol == 1 and coeffs == 0:
                        # p(x) = 1, so p'(x)/p(x) term need not be added
                        presol.append(ybar)
                    # m is a positive integer and there are valid coefficients
                    elif len(coeffs):
                        # Substitute the valid coefficients to get p(x)
                        psol = psol.xreplace(coeffs)
                        # y(x) = ybar(x) + p'(x)/p(x)
                        presol.append(ybar + psol.diff(x)/psol)

    # Remove redundant solutions from the list of existing solutions
    remove = set()
    for i in range(len(presol)):
        for j in range(i+1, len(presol)):
            rem = remove_redundant_sols(presol[i], presol[j], x)
            if rem is not None:
                remove.add(rem)
    sols = [x for x in presol if x not in remove]

    # Step 15 : Inverse transform the solutions of the equation in normal form
    bp = -b2.diff(x)/(2*b2**2) - b1/(2*b2)

    # If general solution is required, compute it from the particular solutions
    if gensol:
        sols = [get_gen_sol_from_part_sol(sols, a, x)]

    # Inverse transform the particular solutions
    presol = [Eq(fx, riccati_inverse_normal(y, x, b1, b2, bp).cancel(extension=True)) for y in sols]
    return presol
