r"""
This module contains :py:meth:`~sympy.solvers.ode.riccati.solve_riccati`,
a function which gives rational solutions to first order Riccati ODEs.
A general first order Riccati ODE is given by

.. math:: y' = b_0(x) + b_1(x)w + b_2(x)w^2

where `b_0, b_1` and `b_2` can be arbitrary functions of `x` with `b_2 \ne 0`.

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

The algorithm implemented here was presented in a Ph.D thesis by
N. Thieu Vo. The entire thesis can be found here -
https://www3.risc.jku.at/publications/download/risc_5387/PhDThesisThieu.pdf

We have only implemented the Rational Riccati solver (Pg 78 in Thesis).
Before we proceed towards the implementation of the algorithm, a
few definitions to understand -

1. Valuation of a Rational Function at `\infty`:
    The valuation of a rational function `p(x)` at `\infty` is equal
    to the difference between the degree of the denominator and the
    numerator of `p(x)`.

    Note: A general definition of valuation of a rational function
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

Step 6 : Find `n` c-vectors (one for each pole) and 1 d-vector using
``construct_c`` and ``construct_d``. Now, determine all the ``2**(n + 1)``
combinations of choosing between 2 choices for each of the `n` c-vectors
and 1 d-vector.

For each of these above combinations, do

Step 8 : Compute `m` in ``compute_degree``

Step 9 : In ``compute_degree``, compute ybar as well. If ybar contains any
term like `\frac{c_{ij}}{(x - oo)^j}` where `c_{ij}` is non-zero, discard
this solution.

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

from sympy.core import S
from sympy.core.numbers import oo
from sympy.core.relational import Eq
from sympy.core.symbol import Symbol, symbols, Dummy
from sympy.functions import sqrt
from sympy.functions.elementary.complexes import sign
from sympy.polys.domains import ZZ
from sympy.polys.polytools import Poly
from sympy.polys.polyroots import roots
from sympy.solvers.solveset import linsolve
from sympy.tensor.indexed import Indexed, IndexedBase

def riccati_normal(w, x, b1, b2):
    # y(x) = -b2(x)*w(x) - b2'(x)/(2*b2(x)) - b1(x)/2
    return -b2*w - b2.diff(x)/(2*b2) - b1/2


def riccati_inverse_normal(y, x, b1, b2, bp=None):
    if bp is None:
        bp = -b2.diff(x)/(2*b2**2) - b1/(2*b2)
    # w(x) = -y(x)/b2(x) - b2'(x)/(2*b2(x)^2) - b1(x)/(2*b2(x))
    return -y/b2 + bp


def linsolve_dict(eq, syms, var):
    sol = linsolve(eq, syms, var)
    if not sol:
        return {}
    return {k:v for k, v in zip(syms, list(sol)[0])}


def find_poles(num, den, x):
    # All roots of denominator are poles of a(x)
    p = roots(den, x)
    # Substitute 1/x for x and check if oo is a pole
    num, den = inverse_transform_poly(num, den, x)

    # Find the poles of a(1/x)
    p_inv = roots(den, x)

    # If 0 is a pole of a(1/x), oo is a pole of a(x)
    if 0 in p_inv:
        p.update({oo: p_inv[0]})
    return p


def val_at_inf(num, den, x):
    # Valuation of a rational function at oo = deg(denom) - deg(numer)
    return den.degree(x) - num.degree(x)


def check_necessary_conds(val_inf, muls):
    return val_inf%2 != 0 or not all([mul == 1 or (mul%2 == 0 and mul >= 2) for mul in muls])


def inverse_transform_poly(num, den, x):
    # Declare for reuse
    one = Poly(1, x)
    xpoly = Poly(x, x)

    # Check if degree of numerator is same as denominator
    pwr = den.degree(x) - num.degree(x)
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
    pwr = num.degree(x) - den.degree(x)
    if pwr > 0:
        return oo*sign(num.LC()/den.LC())
    elif pwr == 0:
        return num.LC()/den.LC()
    else:
        return 0


def construct_c(num, den, x, poles, muls):
    c = []
    for pole, mul in zip(poles, muls):
        c.append([])

        # Case 3
        if mul == 1:
            # If multiplicity is 1, the coefficient to be added
            # in the c-vector is 1 (no choice)
            c[-1].extend([[1], [1]])

        # Case 1
        elif mul == 2:
            # Find the coefficient of 1/(x - pole)**2 in the
            # Laurent series expansion of a(x) about pole.
            if pole != oo:
                sgn, num1, den1 = (num*(x - pole)**2).cancel(den)
                r = (sgn*num1.subs(x, pole))/(den1.subs(x, pole))
            else:
                r = 0

            # If multiplicity is 2, the coefficient to be added
            # in the c-vector is c = (1 +- sqrt(1 + 4*r))/2
            c[-1].extend([[(1 + sqrt(1 + 4*r))/2], [(1 - sqrt(1 + 4*r))/2]])

        # Case 2
        else:
            # Generate the coefficients using the recurrence
            # relation mentioned in (5.14) in the thesis (Pg 80)

            # r_i = mul/2
            ri = mul//2

            # Find the Laurent series coefficients about the pole
            ser = rational_laurent_series(num, den, x, pole, mul, 6)

            # Start with an empty memo to store the coefficients
            # This is for the plus case
            temp = [0 for i in range(ri)]

            # Base Case
            temp[ri-1] = sqrt(ser[0])

            # Iterate backwards to find all coefficients
            for s in range(ri-1, 0, -1):
                sm = 0
                for j in range(s+1, ri):
                    sm += temp[j-1]*temp[ri+s-j-1]
                if s!= 1:
                    temp[s-1] = (ser[mul-ri-s] - sm)/(2*temp[ri-1])

            # Memo for the minus case
            temp1 = list(map(lambda x: -x, temp))

            # Find the 0th coefficient in the recurrence
            temp[0] = (ser[mul-ri-s] - sm + ri*temp[ri-1])/(2*temp[ri-1])
            temp1[0] = (ser[mul-ri-s] - sm + ri*temp1[ri-1])/(2*temp1[ri-1])

            # Add both the plus and minus cases' coefficients
            c[-1].extend([temp, temp1])
    return c


def construct_d(num, den, x, val_inf, mul):
    N = -val_inf//2
    ser = rational_laurent_series(num, den, x, oo, mul, 1)

    # Case 4
    if val_inf < 0:
        temp = [0 for i in range(N+2)]
        temp[N] = sqrt(ser[-mul - 1 + 2*N])
        for s in range(N-1, -2, -1):
            sm = 0
            for j in range(s+1, N):
                sm += temp[j]*temp[N+s-j]
            if s != -1:
                temp[s] = (ser[-mul - 1 + N+s] - sm)/(2*temp[N])
        temp1 = list(map(lambda x: -x, temp))
        temp[-1] = (ser[-mul - 1 + N+s] - temp[N] - sm)/(2*temp[N])
        temp1[-1] = (ser[-mul - 1 + N+s] - temp1[N] - sm)/(2*temp1[N])
        d = [temp, temp1]

    # Case 5
    elif val_inf == 0:
        # List to store coefficients for plus case
        temp = [0, 0]

        # d_0  = sqrt(a_0)
        temp[0] = sqrt(ser[-mul - 1])

        # d_(-1) = a_(-1)/(2*d_0)
        temp[-1] = ser[-mul - 2]/(2*temp[0])

        # Coefficients for the minus case are just the negative
        # of the coefficients for the positive case. Append both.
        d = [temp, list(map(lambda x: -x, temp))]

    # Case 6
    else:
        # s_oo = lim x->0 1/x**2 * a(1/x) which is equivalent to
        # s_oo = lim x->oo x**2 * a(x)
        s_inf = limit_at_inf(Poly(x**2, x)*num, den, x)

        # d_(-1) = (1 +- sqrt(1 + 4*s_oo))/2
        d = [[(1 + sqrt(1 + 4*s_inf))/2], [(1 - sqrt(1 + 4*s_inf))/2]]
    return d


def rational_laurent_series(num, den, x, x0, mul, n):
    m = Symbol('m')
    one = Poly(1, x)
    reverse = False
    if x0 == oo:
        num, den = inverse_transform_poly(num, den, x)
        reverse = True
        x0 = 0
    if x0:
        num = num.transform(Poly(x + x0, x), one)
        den = den.transform(Poly(x + x0, x), one)
        num, den = num.cancel(den, include=True)
    num, den = (num*x**mul).cancel(den, include=True)

    c = IndexedBase('c')
    indices = [c[i] for i in range(max(den.degree(x), num.degree(x) + 1))]
    out, sums = Poly(0, x, domain=ZZ[indices]), 0
    for pw, cf in den.as_dict().items():
        pw = pw[0]
        for i in range(pw, max(den.degree(x), num.degree(x) + 1)):
            out += Poly(cf*c[i - pw]*x**i, x, domain=ZZ[indices])
        sums += cf*c[m - pw]

    if out:
        coeffs = linsolve_dict((num - out).all_coeffs(), list(out.atoms(Indexed)), x)
        for i in range(len(coeffs), n + mul):
            ai = linsolve_dict([sums.subs(m, i).subs(coeffs)], [c[i]], x)
            if type(ai) == dict:
                coeffs.update({c[i]: list(ai.values())[0]})
        coeffs = list(coeffs.items())
        coeffs.sort(key = lambda x: x[0].indices[0])
        coeffs = list(map(lambda x: x[1], coeffs))
    else:
        coeffs = list(map(lambda x: x/den, num.as_dict().values()))
    return coeffs[::-1] if reverse else coeffs


def compute_degree(x, poles, choice, c, d, N):
    ybar = 0
    m = S(d[choice[0]][-1])

    # Calculate the first (nested) summation for ybar
    # as given in Step 9 of the Thesis (Pg 82)
    for i in range(len(poles)):
        for j in range(len(c[i][choice[i + 1]])):
            # If one of the poles is infinity and its coefficient is
            # not zero, the given solution is invalid as there will be
            # a c/(x - oo)^j term in ybar for some c and j
            if poles[i] == oo and c[i][choice[i + 1]][j] != 0:
                return m, ybar, ybar, ybar, False
            ybar += c[i][choice[i + 1]][j]/(x - poles[i])**(j+1)
        m -= c[i][choice[i + 1]][0]

    # Calculate the second summation for ybar
    for i in range(N+1):
        ybar += d[choice[0]][i]*x**i
    numy, deny = [Poly(e, x) for e in ybar.together().as_numer_denom()]
    return m, ybar, numy, deny, True


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
        return psol, linsolve_dict(auxeq.all_coeffs(), psyms, x), True
    else:
        # m == 0 . Check if 1 (x**0) is a solution to the auxiliary equation
        return S(1), auxeq, auxeq == 0


def solve_riccati(eq, fx, x, b0, b1, b2):
    # Step 1 : Convert to Normal Form
    a = -b0*b2 + b1**2/4 - b1.diff(x)/2 + 3*b2.diff(x)**2/(4*b2**2) + b1*b2.diff(x)/(2*b2) - b2.diff(x, 2)/(2*b2)
    a_t = a.together()
    num, den = map(lambda e: Poly(e, x), a_t.as_numer_denom())
    num, den = num.cancel(den, include=True)

    # Step 2
    presol = []

    # Step 3 : a(x) is 0
    if num == 0:
        presol.append(1/(x + get_numbered_constants(eq, 1)))

    # Step 4 : a(x) is a non-zero constant
    elif x not in num.free_symbols.union(den.free_symbols):
        presol.extend([sqrt(a), -sqrt(a)])

    # Step 5 : Find poles and valuation at infinity
    poles = find_poles(num, den, x)
    inf_mul = poles.get(oo, 0)
    poles, muls = list(poles.keys()), list(poles.values())
    val_inf = val_at_inf(num, den, x)

    if len(poles):
        # Check necessary conditions (outlined in the module docstring)
        if check_necessary_conds(val_inf, muls):
            raise ValueError("Rational Solution doesn't exist")

        # Step 6
        # Construct c-vectors for each singular point
        c = construct_c(num, den, x, poles, muls)

        # Construct d vectors for each singular point
        d = construct_d(num, den, x, val_inf, inf_mul)

        # Step 7 : Iterate over all possible combinations and return solutions
        for it in range(2**(len(poles) + 1)):
            # For each possible combination, generate an array of 0's and 1's
            # where 0 means pick 1st choice and 1 means pick the second choice.
            choice = list(map(lambda x: int(x), bin(it)[2:].zfill(len(poles) + 1)))

            # Step 8 and 9 : Compute m and ybar
            m, ybar, numy, deny, exists = compute_degree(x, poles, choice, c, d, -val_inf//2)

            # Step 10 : Check if a valid solution exists. If yes, also check
            # if m is a non-negative integer
            if exists and S(m).is_nonnegative == True and S(m).is_integer == True:

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
                        psol = psol.subs(coeffs)
                        # y(x) = ybar(x) + p'(x)/p(x)
                        presol.append(ybar + psol.diff(x)/psol)

            # If m is an expression with symbols, we generate a Piecewise solution where
            # the solution will exist only if the expression is a nonnegative integer.
            # NOTE: In this case, it seems like x**m + C1 is always a solution
            # for the auxiliary equation. It is however NOT mentioned in the thesis.
            elif len(m.free_symbols):
                C1 = Symbol('C1')
                psol = x**m + C1
                presol.append(ybar + psol.diff(x)/psol)
    # Step 15 : Inverse transform the solutions of the equation in normal form
    bp = -b2.diff(x)/(2*b2**2) - b1/(2*b2)
    sol = [Eq(fx, riccati_inverse_normal(y, x, b1, b2, bp)) for y in presol]
    return sol


# Avoid circular import:
from .ode import get_numbered_constants
