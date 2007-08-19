"""

"""

from sympy.core.basic import Basic, S
from sympy.core.symbol import Symbol
from sympy.core.numbers import Rational
from sympy.core.add import Add
from sympy.core.mul import Mul

from sympy.simplify import simplify, hypersimp, hypersimilar, collect
from sympy.solvers import solve, solve_undetermined_coeffs
from sympy.polynomials import quo, gcd, roots, resultant
from sympy.concrete import nni_roots, product
from sympy.matrices import Matrix, casoratian

def rsolve_poly(coeffs, f, n):
    """Find polynomial solutions to inhomogeneous linear difference
       equation of finite order with polynomial coefficients and
       polynomial inhomogeneous part.

       Formally, given a linear difference equation Ly = f where

          L = p_r * E^k + ... + p_1 * E + p_0  ;  E^k y = y(n+k)

       with p_i : i=0..r and 'f' being polynomialsin 'n' we seek for
       solutions y(n) in K[n] where 'K' is a field of characteristic
       zero. In this implementation we have K = C.

       The algorithm performs two basic steps:

           (1) Compute degree N of the general polynomial solution.
           (2) Find all polynomials of degree N or less of Ly = f.

       There are two methods for computing the polynomial solutions.
       If the degree bound is relatively small, ie. it is smaller than
       or equal to the order of the recurrence, then naive method of
       undetermined coefficients is being used. This gives system
       of algebraic equations with N+1 unknowns.

       In the other case, the algorithm performs transformation of the
       initial equation to an equivalent one, for which the system of
       algebraic equations has only 'r' undeterminates. This method is
       quite sophisticated (in comparison with the naive one) and was
       invented together by Abramov, Bronstein and Petkovsek.

       It is possible to generalize the algorithm implemented here to
       the case of linear q-difference and differential equations.

       Note the fact that this procedure is highly specific and does
       not perform any analysis of the arguments given. If not sure
       use rsolve() instead, which will run appropriate algorithm
       for the given problem.

       If there is need to run this procedure directly then the list of
       coefficients must be in the following form: [ p_0, p_1, ..., p_r ].

       Lets say that we would like to compute m-th Bernoulli polynomial
       up to a constant. For this we can use b(n+1) - b(n) == m*n**(m-1)
       recurrence, which has solution b(n) = B_m + C. For example:

       >>> from sympy.core import Symbol
       >>> n = Symbol('n', integer=True)

       >>> rsolve_poly([-1, 1], 4*n**3, n)
       C0 + n**2 + n**4 - 2*n**3

       For more information on the implemented algorithm refer to:

       [1] S. A. Abramov, M. Bronstein and M. Petkovsek, On polynomial
           solutions of linear operator equations, in: T. Levelt, ed.,
           Proc. ISSAC '95, ACM Press, New York, 1995, 290-296.

       [2] M. Petkovsek, Hypergeometric solutions of linear recurrences
           with polynomial coefficients, J. Symbolic Computation,
           14 (1992), 243-264.

       [3] M. Petkovsek, H. S. Wilf, D. Zeilberger, A = B, 1996.

    """
    f = Basic.sympify(f)

    if not f.is_polynomial(n):
        return None

    homogeneous = isinstance(f, Basic.Zero)

    if not homogeneous:
        df = f.as_polynomial(n).degree()

    coeffs = map(Basic.sympify, coeffs)

    r = len(coeffs)-1

    polys = [Basic.Zero()]*(r+1)
    terms = [ (Basic.Zero(), Basic.NegativeInfinity()) ]*(r+1)

    for i in xrange(0, r+1):
        for j in xrange(i, r+1):
            polys[i] += S.Binomial(j, i)*coeffs[j]

        polys[i] = polys[i].expand()

        if not isinstance(polys[i], Basic.Zero):
            terms[i] = polys[i].as_polynomial(n).coeffs[0]

    d = b = terms[0][1]

    for i in xrange(1, r+1):
        if terms[i][1] > d:
            d = terms[i][1]

        if terms[i][1] - i > b:
            b = terms[i][1] - i

    d, b = int(d), int(b)

    x = Symbol('x', dummy=True)

    degree_poly = Basic.Zero()

    for i in xrange(0, r+1):
        if terms[i][1] - i == b:
            degree_poly += terms[i][0]*S.FallingFactorial(x, i)

    _nni_roots = nni_roots(degree_poly, x)

    if _nni_roots != []:
        N = [max(_nni_roots)]
    else:
        N = []

    if homogeneous:
        N += [-b-1]
    else:
        N += [df-b, -b-1]

    N = int(max(N))

    if N < 0:
        return None

    if N <= r:
        C = []
        y = E = Basic.Zero()

        for i in xrange(0, N+1):
            C.append(Symbol('C'+str(i)))
            y += C[i] * n**i

        for i in xrange(0, r+1):
            E += coeffs[i]*y.subs(n, n+i)

        return y.subs_dict(solve_undetermined_coeffs(E-f, C, n))
    else:
        A = r
        U = N+A+b+1

        _nni_roots = nni_roots(polys[r], n)

        if _nni_roots != []:
            a = max(_nni_roots) + 1
        else:
            a = Basic.Zero()

        def zero_vector(k):
            return [Basic.Zero()] * k

        def one_vector(k):
            return [Basic.One()] * k

        def delta(p, k):
            B = Basic.One()
            D = p.subs(n, a+k)

            for i in xrange(1, k+1):
                B *= -Rational(k-i+1, i)
                D += B * p.subs(n, a+k-i)

            return D

        alpha = {}

        for i in xrange(-A, d+1):
            I = one_vector(d+1)

            for k in xrange(1, d+1):
                I[k] = I[k-1] * (x+i-k+1)/k

            alpha[i] = Basic.Zero()

            for j in xrange(0, A+1):
                for k in xrange(0, d+1):
                    B = S.Binomial(k, i+j)
                    D = delta(polys[j], k)

                    alpha[i] += I[k]*B*D

        V = Matrix(U, A, lambda i, j: int(i == j))

        if homogeneous:
            for i in xrange(A, U):
                v = zero_vector(A)

                for k in xrange(1, A+b+1):
                    if i - k < 0:
                        break

                    B = alpha[k-A].subs(x, i-k)

                    for j in xrange(0, A):
                        v[j] += B * V[i-k, j]

                denom = alpha[-A].subs(x, i)

                for j in xrange(0, A):
                    V[i, j] = -v[j] / denom
        else:
            G = zero_vector(U)

            for i in xrange(A, U):
                v = zero_vector(A)
                g = Basic.Zero()

                for k in xrange(1, A+b+1):
                    if i - k < 0:
                        break

                    B = alpha[k-A].subs(x, i-k)

                    for j in xrange(0, A):
                        v[j] += B * V[i-k, j]

                    g += B * G[i-k]

                denom = alpha[-A].subs(x, i)

                for j in xrange(0, A):
                    V[i, j] = -v[j] / denom

                G[i] = (delta(f, i-A) - g) / denom

        P, Q = one_vector(U), zero_vector(A)

        for i in xrange(1, U):
            P[i] = (P[i-1] * (n-a-i+1)/i).expand()

        for i in xrange(0, A):
            Q[i] = Add(*[ (v*p).expand() for v, p in zip(V[:,i], P) ])

        if not homogeneous:
            h = Add(*[ (g*p).expand() for g, p in zip(G, P) ])

        C = [ Symbol('C'+str(i)) for i in xrange(0, A) ]

        g = lambda i: Add(*[ c*delta(q, i) for c, q in zip(C, Q) ])

        if homogeneous:
            E = [ g(i) for i in xrange(N+1, U) ]
        else:
            E = [ g(i) + delta(h, i) for i in xrange(N+1, U) ]

        if E != []:
            solutions = solve(E, C)
        else:
            solutions = {}

        result = []

        for c, s in zip(C, Q):
            if c in solutions:
                result.append((solutions[c] * q).expand())
            else:
                result.append((c * s).expand())

        if homogeneous:
            return Add(*result)
        else:
            return Add(*result) + h

def rsolve_ratio(coeffs, f, n):
    """Find rational solutions to inhomogeneous linear difference
       equation of finite order with polynomial coefficients and
       polynomial inhomogeneous part.

       Formally, given a linear difference equation Ly = f where

          L = p_r * E^k + ... + p_1 * E + p_0  ;  E^k y = y(n+k)

       with p_i : i=0..r and 'f' being polynomialsin 'n', we seek
       for solutions y(n) of the form u(n)/v(n). Of course both
       u and v are polynomials over a field of characteristic
       zero. In this implementation we have K = C.

       The algorithm performs two basic steps:

           (1) Compute polynomial v(n) which can be used as universal
               denominator of any rational solution of equation Ly = f.

           (2) Construct new linear difference equation by substitution
               y(n) = u(n)/v(n) and solve it for u(n) finding all its
               polynomial solutions. Return None if none were found.

       Algorithm implemented here is a revised version of the original
       Abramov's algorithm, developed in 1989. The new approach is much
       simpler to implement and has better overall efficiency. This
       method can be easily adapted to q-difference equations case.

       For more information on the implemented algorithm refer to:

       [1] S. A. Abramov, Rational solutions of linear difference
           and q-difference equations with polynomial coefficients,
           in: T. Levelt, ed., Proc. ISSAC '95, ACM Press, New York,
           1995, 285-289

    """
    f = Basic.sympify(f)

    if not f.is_polynomial(n):
        return None

    coeffs = map(Basic.sympify, coeffs)

    r = len(coeffs)-1

    A, B = coeffs[r], coeffs[0]
    A = A.subs(n, n-r).expand()

    h = Symbol('h', dummy=True)

    res = resultant(A, B.subs(n, n+h), n)

    if not res.is_polynomial(h):
        p, q = res.as_numer_denom()
        res = quo(p, q, h)

    _nni_roots = nni_roots(res, h)

    if _nni_roots == []:
        return rsolve_poly(coeffs, f, n, **hints)
    else:
        C, numers = S.One, [S.Zero]*(r+1)

        for i in xrange(int(max(_nni_roots)), -1, -1):
            d = gcd(A, B.subs(n, n+i), n)

            A = quo(A, d, n)
            B = quo(B, d.subs(n, n-i), n)

            C *= Mul(*[ d.subs(n, n-j) for j in xrange(0, i+1) ])

        denoms = [ C.subs(n, n+i) for i in range(0, r+1) ]

        for i in range(0, r+1):
            g = gcd(coeffs[i], denoms[i], n)

            numers[i] = quo(coeffs[i], g, n)
            denoms[i] = quo(denoms[i], g, n)

        for i in xrange(0, r+1):
            numers[i] *= Mul(*(denoms[:i] + denoms[i+1:]))

        result = rsolve_poly(numers, f * Mul(*denoms), n, **hints)

        if result is not None:
            if hints.get('symbols', False):
                return (simplify(result[0] / C), result[1])
            else:
                return simplify(result / C)
        else:
            return None

def rsolve_hyper(coeffs, f, n, **hints):
    """Find hypergeometric solutions to inhomogeneous linear
       difference equation of finite order with polynomial
       coefficients.

       Formally, given a linear difference equation Ly = f where

         L = p_r * E^k + ... + p_1 * E + p_0  ;  E^k y = y(n+k)

       with p_i : i=0..r being polynomials in 'n' and where the
       inhomogeneous part is hypergeometric or is a combination
       of pairwise dissimilar hypergeometric terms, we seek for
       all hypergeometric.

       Hypergeometric terms are integer sequences anihilated by
       first order linear difference equations with polynomial
       coefficients or, in simpler words, if consecutive term
       ratio is quotient of polynomials. This gives us notion
       of closed form.

       Class of D'Alembertian terms is much richer and contains
       also sequences which are products of hypergeometric and
       not definitely summable. ie. not Gosper summable terms.
       Classical example of usage of this generalisation is
       formula for number of derangements (permutations in
       which none of the elements is in proper position).

       For efficiency reason all computations are done in field
       of real numbers. However if complex solutions of equations
       with constant coefficients are needed then use method of
       characteristic polynomials instead.

       For more information on the implemented algorithm refer to:

       [1] M. Petkovsek, Hypergeometric solutions of linear recurrences
           with polynomial coefficients, J. Symbolic Computation,
           14 (1992), 243-264.

       [2] M. Petkovsek, H. S. Wilf, D. Zeilberger, A = B, 1996.

    """
    coeffs = map(Basic.sympify, coeffs)

    f = Basic.sympify(f)

    r, kernel = len(coeffs)-1, []

    if not isinstance(f, Basic.Zero):
        if f.is_hypergeometric(n):
            inhomogeneous = [f]
        else:
            if isinstance(f, Basic.Add):
                similar = {}

                for g in f.expand():
                    if not g.is_hypergeometric(n):
                        return None

                    for h in similar.iterkeys():
                        if hypersimilar(g, h, n):
                            similar[h] += g
                            break
                    else:
                        similar[g] = S.Zero

                inhomogeneous = []

                for g, h in similar.iteritems():
                    inhomogeneous.append(g+h)
            else:
                return None

        for i, g in enumerate(inhomogeneous):
            coeff, polys = S.One, coeffs[:]

            s = hypersimp(g, n)

            for j in range(1, r+1):
                coeff *= s.subs(n, n+j-1)
                polys[j] *= coeff

            R = rsolve_ratio(polys, S.One, n)

            if not (R is None or isinstance(R, Basic.Zero)):
                inhomogeneous[i] *= R
            else:
                return None

            result = Add(*inhomogeneous)
    else:
        result = S.Zero

    Z = Symbol('Z', dummy=True)

    p, q = coeffs[0], coeffs[r].subs(n, n-r+1)

    p_factors = [ z for z in set(roots(p, n)) if z.is_real ]
    q_factors = [ z for z in set(roots(q, n)) if z.is_real ]

    factors = [ (S.One, S.One) ]

    for p in p_factors:
        for q in q_factors:
            if p.is_integer and q.is_integer and p <= q:
                continue
            else:
                factors += [(n-p, n-q)]

    p = [ (n-p, S.One) for p in p_factors ]
    q = [ (S.One, n-q) for q in q_factors ]

    factors = p + factors + q

    for A, B in factors:
        polys, degrees = [], []
        D = A*B.subs(n, n+r-1)

        for i in range(0, r+1):
            a = Mul(*[ A.subs(n, n+j) for j in range(0, i) ])
            b = Mul(*[ B.subs(n, n+j) for j in range(i, r) ])

            poly = (quo(coeffs[i]*a*b, D, n))
            polys.append(poly.as_polynomial(n))

            if not isinstance(poly, Basic.Zero):
                degrees.append(polys[i].degree())

        d, poly = max(degrees), S.Zero

        for i in range(0, r+1):
            coeff = polys[i].nth_coeff(d)

            if not isinstance(coeff, Basic.Zero):
                poly += coeff * Z**i

        for z in set(roots(poly, Z)):
            if not z.is_real or z.is_zero:
                continue

            C = rsolve_poly([ polys[i]*z**i for i in range(r+1) ], 0, n)

            if C is not None and not isinstance(C, Basic.Zero):
                K = simplify(z*(A*C.subs(n, n+1))/(B*C))

                if casoratian(kernel+[K], n) != 0:
                    kernel.append(K)

    symbols = [ Symbol('C'+str(i)) for i in xrange(len(kernel)) ]

    for C, ker in zip(symbols, kernel):
        result += C * product(ker, (n, 0, n-1))

    if hints.get('symbols', False):
        return (result, symbols)
    else:
        return result

def rsolve(eq, seq):
    """

    """
    pass
