"""
Expand Hypergeometric (and Meijer G) functions into named
special functions.

The algorithm for doing this uses a collection of lookup tables of
hypergeometric functions, and various of their properties, to expand
many hypergeometric functions in terms of special functions.

It is based on the following paper:
      Kelly B. Roach.  Meijer G Function Representations.
      In: Proceedings of the 1997 International Symposium on Symbolic and
      Algebraic Computation, pages 205-211, New York, 1997. ACM.

It is described in great(er) detail in the Sphinx documentation.
"""
# SUMMARY OF EXTENSIONS FOR MEIJER G FUNCTIONS
#
# o z**rho G(ap, bq; z) = G(ap + rho, bq + rho; z)
#
# o denote z*d/dz by D
#
# o It is helpful to keep in mind that ap and bq play essentially symmetric
#   roles: G(1/z) has slightly altered parameters, with ap and bq interchanged.
#
# o There are four shift operators:
#   A_J = b_J - D,     J = 1, ..., n
#   B_J = 1 - a_j + D, J = 1, ..., m
#   C_J = -b_J + D,    J = m+1, ..., q
#   D_J = a_J - 1 - D, J = n+1, ..., p
#
#   A_J, C_J increment b_J
#   B_J, D_J decrement a_J
#
# o The corresponding four inverse-shift operators are defined if there
#   is no cancellation. Thus e.g. an index a_J (upper or lower) can be
#   incremented if a_J != b_i for i = 1, ..., q.
#
# o Order reduction: if b_j - a_i is a non-negative integer, where
#   j <= m and i > n, the corresponding quotient of gamma functions reduces
#   to a polynomial. Hence the G function can be expressed using a G-function
#   of lower order.
#   Similarly if j > m and i <= n.
#
#   Secondly, there are paired index theorems [Adamchik, The evaluation of
#   integrals of Bessel functions via G-function identities]. Suppose there
#   are three parameters a, b, c, where a is an a_i, i <= n, b is a b_j,
#   j <= m and c is a denominator parameter (i.e. a_i, i > n or b_j, j > m).
#   Suppose further all three differ by integers.
#   Then the order can be reduced.
#   TODO work this out in detail.
#
# o An index quadruple is called suitable if its order cannot be reduced.
#   If there exists a sequence of shift operators transforming one index
#   quadruple into another, we say one is reachable from the other.
#
# o Deciding if one index quadruple is reachable from another is tricky. For
#   this reason, we use hand-built routines to match and instantiate formulas.
#
from sympy.core import S, Dummy, symbols, sympify, Tuple, expand, I, Mul
from sympy.core.mod import Mod
from sympy.functions.special.hyper import hyper
from sympy.utilities.misc import default_sort_key
from sympy import SYMPY_DEBUG

from sympy.utilities.timeutils import timethis
_timeit = timethis('meijerg')

# leave add formulae at the top for easy reference


def add_formulae(formulae):
    """ Create our knowledge base. """

    a, b, c, z = symbols('a b c, z', cls=Dummy)

    def add(ap, bq, res):
        formulae.append(Formula(ap, bq, z, res, (a, b, c)))

    def addb(ap, bq, B, C, M):
        formulae.append(Formula(ap, bq, z, None, (a, b, c), B, C, M))

    # Luke, Y. L. (1969), The Special Functions and Their Approximations,
    # Volume 1, section 6.2

    from sympy import (exp, sqrt, root, cosh, log, asin, atan, I, lowergamma, cos,
                       atanh, besseli, gamma, erf, pi, sin, besselj, Ei,
                       EulerGamma, Shi, sinh, cosh, Chi, diag, Matrix,
                       fresnels, fresnelc)
    from sympy.functions.special.hyper import (HyperRep_atanh,
        HyperRep_power1, HyperRep_power2, HyperRep_log1, HyperRep_asin1,
        HyperRep_asin2, HyperRep_sqrts1, HyperRep_sqrts2, HyperRep_log2,
        HyperRep_cosasin, HyperRep_sinasin)
    from sympy import polar_lift, exp_polar

    # 0F0
    add((), (), exp(z))

    # 1F0
    add((-a, ), (), HyperRep_power1(a, z))

    # 2F1
    addb((a, a - S.Half), (2*a, ),
         Matrix([HyperRep_power2(a, z),
                 HyperRep_power2(a + S(1)/2, z)/2]),
         Matrix([[1, 0]]),
         Matrix([[(a-S.Half)*z/(1 - z), (S.Half - a)*z/(1 - z)],
                 [a/(1 - z), a*(z - 2)/(1 - z)]]))
    addb((1, 1), (2, ),
         Matrix([HyperRep_log1(z), 1]), Matrix([[-1/z, 0]]),
         Matrix([[0, z/(z - 1)], [0, 0]]))
    addb((S.Half, 1), (S('3/2'), ),
         Matrix([HyperRep_atanh(z), 1]),
         Matrix([[1, 0]]),
         Matrix([[-S(1)/2, 1/(1 - z)/2], [0, 0]]))
    addb((S.Half, S.Half), (S('3/2'), ),
         Matrix([HyperRep_asin1(z), HyperRep_power1(-S(1)/2, z)]),
         Matrix([[1, 0]]),
         Matrix([[-S(1)/2, S(1)/2], [0, z/(1 - z)/2]]))
    addb((-a, S.Half - a), (S.Half, ),
         Matrix([HyperRep_sqrts1(a, z), -HyperRep_sqrts2(a - S(1)/2, z)]),
         Matrix([[1, 0]]),
         Matrix([[0, a],
                 [z*(2*a - 1)/2/(1 - z), S.Half - z*(2*a - 1)/(1 - z)]]))

    # A. P. Prudnikov, Yu. A. Brychkov and O. I. Marichev (1990).
    # Integrals and Series: More Special Functions, Vol. 3,.
    # Gordon and Breach Science Publisher
    addb([a, -a], [S.Half],
         Matrix([HyperRep_cosasin(a, z), HyperRep_sinasin(a, z)]),
         Matrix([[1, 0]]),
         Matrix([[0, -a], [a*z/(1 - z), 1/(1 - z)/2]]))
    addb([1, 1], [3*S.Half],
         Matrix([HyperRep_asin2(z), 1]), Matrix([[1, 0]]),
         Matrix([[(z - S.Half)/(1 - z), 1/(1 - z)/2], [0, 0]]))

    # 3F2
    addb([-S.Half, 1, 1], [S.Half, 2],
         Matrix([z*HyperRep_atanh(z), HyperRep_log1(z), 1]),
         Matrix([[-S(2)/3, -S(1)/(3*z), S(2)/3]]),
         Matrix([[S(1)/2, 0, z/(1 - z)/2],
                 [0, 0, z/(z - 1)],
                 [0, 0, 0]]))
    # actually the formula for 3/2 is much nicer ...
    addb([-S.Half, 1, 1], [2, 2],
         Matrix([HyperRep_power1(S(1)/2, z), HyperRep_log2(z), 1]),
         Matrix([[S(4)/9 - 16/(9*z), 4/(3*z), 16/(9*z)]]),
         Matrix([[z/2/(z - 1), 0, 0], [1/(2*(z - 1)), 0, S.Half], [0, 0, 0]]))

    # 1F1
    addb([1], [b], Matrix([z**(1 - b) * exp(z) * lowergamma(b - 1, z), 1]),
         Matrix([[b - 1, 0]]), Matrix([[1 - b + z, 1], [0, 0]]))
    addb([a], [2*a],
         Matrix([z**(S.Half - a)*exp(z/2)*besseli(a - S.Half, z/2)
                 * gamma(a + S.Half)/4**(S.Half - a),
                 z**(S.Half - a)*exp(z/2)*besseli(a + S.Half, z/2)
                 * gamma(a + S.Half)/4**(S.Half - a)]),
         Matrix([[1, 0]]),
         Matrix([[z/2, z/2], [z/2, (z/2 - 2*a)]]))
    mz = polar_lift(-1)*z
    addb([a], [a + 1],
         Matrix([mz**(-a)*a*lowergamma(a, mz), a*exp(z)]),
         Matrix([[1, 0]]),
         Matrix([[-a, 1], [0, z]]))
    # This one is redundant.
    add([-S.Half], [S.Half], exp(z) - sqrt(pi*z)*(-I)*erf(I*sqrt(z)))

    # Added to get nice results for Laplace transform of Fresnel functions
    # http://functions.wolfram.com/07.22.03.6437.01
    # Basic rule
    #add([1], [S(3)/4, S(5)/4],
    #    sqrt(pi) * (cos(2*sqrt(polar_lift(-1)*z))*fresnelc(2*root(polar_lift(-1)*z,4)/sqrt(pi)) +
    #                sin(2*sqrt(polar_lift(-1)*z))*fresnels(2*root(polar_lift(-1)*z,4)/sqrt(pi)))
    #    / (2*root(polar_lift(-1)*z,4)))
    # Manually tuned rule
    addb([1], [S(3)/4, S(5)/4],
         Matrix([ sqrt(pi)*(I*sinh(2*sqrt(z))*fresnels(2*root(z,4)*exp(I*pi/4)/sqrt(pi))
                            + cosh(2*sqrt(z))*fresnelc(2*root(z,4)*exp(I*pi/4)/sqrt(pi)))
                  * exp(-I*pi/4)/(2*root(z,4)),
                  sqrt(pi)*root(z,4)*(sinh(2*sqrt(z))*fresnelc(2*root(z,4)*exp(I*pi/4)/sqrt(pi))
                                      + I*cosh(2*sqrt(z))*fresnels(2*root(z,4)*exp(I*pi/4)/sqrt(pi)))
                  *exp(-I*pi/4)/2,
                  1 ]),
         Matrix([[1, 0, 0]]),
         Matrix([[-S(1)/4, 1,      S(1)/4],
                 [ z,      S(1)/4, 0],
                 [ 0,      0,      0]]))

    # 2F2
    addb([S.Half, a], [S(3)/2, a + 1],
         Matrix([a/(2*a - 1)*(-I)*sqrt(pi/z)*erf(I*sqrt(z)),
                 a/(2*a - 1)*(polar_lift(-1)*z)**(-a)*
                 lowergamma(a, polar_lift(-1)*z),
                 a/(2*a - 1)*exp(z)]),
         Matrix([[1, -1, 0]]),
         Matrix([[-S.Half, 0, 1], [0, -a, 1], [0, 0, z]]))
    # We make a "basis" of four functions instead of three, and give EulerGamma
    # an extra slot (it could just be a coefficient to 1). The advantage is
    # that this way Polys will not see multivariate polynomials (it treats
    # EulerGamma as an indeterminate), which is *way* faster.
    addb([1, 1], [2, 2],
         Matrix([Ei(z) - log(z), exp(z), 1, EulerGamma]),
         Matrix([[1/z, 0, 0, -1/z]]),
         Matrix([[0, 1, -1, 0], [0, z, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]))

    # 0F1
    add((), (S.Half, ), cosh(2*sqrt(z)))
    addb([], [b],
         Matrix([gamma(b)*z**((1 - b)/2)*besseli(b - 1, 2*sqrt(z)),
                 gamma(b)*z**(1 - b/2)*besseli(b, 2*sqrt(z))]),
         Matrix([[1, 0]]), Matrix([[0, 1], [z, (1 - b)]]))

    # 0F3
    x = 4*z**(S(1)/4)

    def fp(a, z):
        return besseli(a, x) + besselj(a, x)

    def fm(a, z):
        return besseli(a, x) - besselj(a, x)

    # TODO branching
    addb([], [S.Half, a, a + S.Half],
         Matrix([fp(2*a - 1, z), fm(2*a, z)*z**(S(1)/4),
                 fm(2*a - 1, z)*sqrt(z), fp(2*a, z)*z**(S(3)/4)])
           * 2**(-2*a)*gamma(2*a)*z**((1 - 2*a)/4),
         Matrix([[1, 0, 0, 0]]),
         Matrix([[0, 1, 0, 0],
                 [0, S(1)/2 - a, 1, 0],
                 [0, 0, S(1)/2, 1],
                 [z, 0, 0, 1 - a]]))
    x = 2*(4*z)**(S(1)/4)*exp_polar(I*pi/4)
    addb([], [a, a + S.Half, 2*a],
         (2*sqrt(polar_lift(-1)*z))**(1 - 2*a)*gamma(2*a)**2 *
         Matrix([besselj(2*a - 1, x)*besseli(2*a - 1, x),
                 x*(besseli(2*a, x)*besselj(2*a - 1, x)
                    - besseli(2*a - 1, x)*besselj(2*a, x)),
                 x**2*besseli(2*a, x)*besselj(2*a, x),
                 x**3*(besseli(2*a, x)*besselj(2*a - 1, x)
                       + besseli(2*a - 1, x)*besselj(2*a, x))]),
         Matrix([[1, 0, 0, 0]]),
         Matrix([[0, S(1)/4, 0, 0],
                 [0, (1 - 2*a)/2, -S(1)/2, 0],
                 [0, 0, 1 - 2*a, S(1)/4],
                 [-32*z, 0, 0, 1 - a]]))

    # 1F2
    addb([a], [a - S.Half, 2*a],
         Matrix([z**(S.Half - a)*besseli(a - S.Half, sqrt(z))**2,
                 z**(1 - a)*besseli(a - S.Half, sqrt(z))
                         *besseli(a - S(3)/2, sqrt(z)),
                 z**(S(3)/2 - a)*besseli(a-S(3)/2, sqrt(z))**2]),
         Matrix([[-gamma(a + S.Half)**2/4**(S.Half - a),
                 2*gamma(a - S.Half)*gamma(a + S.Half)/4**(1 - a),
                 0]]),
         Matrix([[1 - 2*a, 1, 0], [z/2, S.Half - a, S.Half], [0, z, 0]]))
    addb([S.Half], [b, 2 - b],
         pi*(1 - b)/sin(pi*b)*
         Matrix([besseli(1 - b, sqrt(z))*besseli(b - 1, sqrt(z)),
                 sqrt(z)*(besseli(-b, sqrt(z))*besseli(b - 1, sqrt(z))
                          + besseli(1 - b, sqrt(z))*besseli(b, sqrt(z))),
                 besseli(-b, sqrt(z))*besseli(b, sqrt(z))]),
         Matrix([[1, 0, 0]]),
         Matrix([[b-1, S(1)/2, 0],
                 [z, 0, z],
                 [0, S(1)/2, -b]]))
    addb([S(1)/2], [S(3)/2, S(3)/2],
         Matrix([Shi(2*sqrt(z))/2/sqrt(z), sinh(2*sqrt(z))/2/sqrt(z),
                 cosh(2*sqrt(z))]),
         Matrix([[1, 0, 0]]),
         Matrix([[-S.Half, S.Half, 0], [0, -S.Half, S.Half], [0, 2*z, 0]]))

    # FresnelS
    # Basic rule
    #add([S(3)/4], [S(3)/2,S(7)/4], 6*fresnels( exp(pi*I/4)*root(z,4)*2/sqrt(pi) ) / ( pi * (exp(pi*I/4)*root(z,4)*2/sqrt(pi))**3 ) )
    # Manually tuned rule
    addb([S(3)/4], [S(3)/2,S(7)/4],
         Matrix([ fresnels( exp(pi*I/4)*root(z,4)*2/sqrt(pi) ) / ( pi * (exp(pi*I/4)*root(z,4)*2/sqrt(pi))**3 ),
                  sinh(2*sqrt(z))/sqrt(z),
                  cosh(2*sqrt(z)) ]),
         Matrix([[6, 0, 0]]),
         Matrix([[-S(3)/4,  S(1)/16, 0],
                 [ 0,      -S(1)/2,  1],
                 [ 0,       z,       0]]))

    # FresnelC
    # Basic rule
    #add([S(1)/4], [S(1)/2,S(5)/4], fresnelc( exp(pi*I/4)*root(z,4)*2/sqrt(pi) ) / ( exp(pi*I/4)*root(z,4)*2/sqrt(pi) ) )
    # Manually tuned rule
    addb([S(1)/4], [S(1)/2,S(5)/4],
         Matrix([ sqrt(pi)*exp(-I*pi/4)*fresnelc(2*root(z,4)*exp(I*pi/4)/sqrt(pi))/(2*root(z,4)),
                  cosh(2*sqrt(z)),
                  sinh(2*sqrt(z))*sqrt(z) ]),
         Matrix([[1, 0, 0]]),
         Matrix([[-S(1)/4,  S(1)/4, 0     ],
                 [ 0,       0,      1     ],
                 [ 0,       z,      S(1)/2]]))

    # 2F3
    # XXX with this five-parameter formula is pretty slow with the current
    #     Formula.find_instantiations (creates 2!*3!*3**(2+3) ~ 3000
    #     instantiations ... But it's not too bad.
    addb([a, a + S.Half], [2*a, b, 2*a - b + 1],
         gamma(b)*gamma(2*a - b + 1) * (sqrt(z)/2)**(1 - 2*a) *
         Matrix([besseli(b - 1, sqrt(z))*besseli(2*a - b, sqrt(z)),
                 sqrt(z)*besseli(b, sqrt(z))*besseli(2*a - b, sqrt(z)),
                 sqrt(z)*besseli(b - 1, sqrt(z))*besseli(2*a - b + 1, sqrt(z)),
                 besseli(b, sqrt(z))*besseli(2*a - b + 1, sqrt(z))]),
         Matrix([[1, 0, 0, 0]]),
         Matrix([[0, S(1)/2, S(1)/2, 0],
                 [z/2, 1 - b, 0, z/2],
                 [z/2, 0, b - 2*a, z/2],
                 [0, S(1)/2, S(1)/2, -2*a]]))
    # (C/f above comment about eulergamma in the basis).
    addb([1, 1], [2, 2, S(3)/2],
         Matrix([Chi(2*sqrt(z)) - log(2*sqrt(z)),
                 cosh(2*sqrt(z)), sqrt(z)*sinh(2*sqrt(z)), 1, EulerGamma]),
         Matrix([[1/z, 0, 0, 0, -1/z]]),
         Matrix([[0, S(1)/2, 0, -S(1)/2, 0],
                 [0, 0, 1, 0, 0],
                 [0, z, S(1)/2, 0, 0],
                 [0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0]]))


def add_meijerg_formulae(formulae):
    from sympy import Matrix, gamma, uppergamma, exp, Si, Ci, sin, cos, sqrt, pi

    a, b, c, z = map(Dummy, 'abcz')
    rho = Dummy('rho')

    def add(an, ap, bm, bq, B, C, M, matcher):
        formulae.append(MeijerFormula(an, ap, bm, bq, z, [a, b, c, rho],
                                      B, C, M, matcher))

    def detect_uppergamma(iq):
        x = iq.an[0]
        y, z = iq.bm
        swapped = False
        if not Mod((x - y).simplify(), 1):
            swapped = True
            (y, z) = (z, y)
        if Mod((x - z).simplify(), 1) or x > z:
            return None
        l = [y, x]
        if swapped:
            l = [x, y]
        return {rho: y, a: x - y}, IndexQuadruple([x], [], l, [])

    add([a + rho], [], [rho, a + rho], [],
        Matrix([gamma(1 - a)*z**rho*exp(z)*uppergamma(a, z),
                gamma(1 - a)*z**(a + rho)]),
        Matrix([[1, 0]]),
        Matrix([[rho+z, -1], [0, a+rho]]),
        detect_uppergamma)

    def detect_3113(iq):
        """http://functions.wolfram.com/07.34.03.0984.01"""
        x = iq.an[0]
        u, v, w = iq.bm
        if Mod((u - v).simplify(), 1) == 0:
            if Mod((v - w).simplify(), 1) == 0:
                return
            sig = (S(1)/2, S(1)/2, S(0))
            x1, x2, y = u, v, w
        else:
            if Mod((x - u).simplify(), 1) == 0:
                sig = (S(1)/2, S(0), S(1)/2)
                x1, y, x2 = u, v, w
            else:
                sig = (S(0), S(1)/2, S(1)/2)
                y, x1, x2 = u, v, w

        if (Mod((x - x1).simplify(), 1) != 0 or
            Mod((x - x2).simplify(), 1) != 0 or
            Mod((x - y).simplify(), 1) != S(1)/2 or
            x > x1 or x > x2):
            return

        return {a: x}, IndexQuadruple([x], [], [x - S(1)/2 + t for t in sig], [])

    s = sin(2*sqrt(z))
    c_ = cos(2*sqrt(z))
    S_ = Si(2*sqrt(z)) - pi/2
    C = Ci(2*sqrt(z))
    add([a], [], [a, a, a - S(1)/2], [],
        Matrix([sqrt(pi)*z**(a - S(1)/2)*(c_*S_ - s*C),
                sqrt(pi)*z**a*(s*S_ + c_*C),
                sqrt(pi)*z**a]),
        Matrix([[-2, 0, 0]]),
        Matrix([[a - S(1)/2, -1, 0], [z, a, S(1)/2], [0, 0, a]]),
        detect_3113)


def make_simp(z):
    """ Create a function that simplifies rational functions in ``z``. """

    def simp(expr):
        """ Efficiently simplify the rational function ``expr``. """
        from sympy import poly
        numer, denom = expr.as_numer_denom()
        c, numer, denom = poly(numer, z).cancel(poly(denom, z))
        return c * numer.as_expr() / denom.as_expr()

    return simp


def debug(*args):
    if SYMPY_DEBUG:
        for a in args:
            print a,
        print


class IndexPair(object):
    """ Holds a pair of indices, and methods to compute their invariants. """

    def __init__(self, ap, bq):
        from sympy import expand, Tuple
        self.ap = Tuple(*map(expand, ap))
        self.bq = Tuple(*map(expand, bq))

    @property
    def sizes(self):
        return (len(self.ap), len(self.bq))

    def __repr__(self):
        return 'IndexPair(%s, %s)' % (self.ap, self.bq)

    def compute_buckets(self, oabuckets=None, obbuckets=None):
        """
        Partition ``ap`` and ``bq`` mod 1.

        Partition parameters ``ap``, ``bq`` into buckets, i.e., return two
        dicts abuckets, bbuckets such that every key in (ab)buckets is a
        rational in the range [0, 1) [represented either by a Rational or a Mod
        object], and such that (ab)buckets[key] is a tuple of all elements of
        respectively ``ap`` or ``bq`` that are congruent to ``key`` modulo 1.

        If oabuckets, obbuckets is specified, try to use the same Mod objects
        for parameters where possible.

        >>> from sympy.simplify.hyperexpand import IndexPair
        >>> from sympy import S
        >>> ap = (S(1)/2, S(1)/3, S(-1)/2, -2)
        >>> bq = (1, 2)
        >>> IndexPair(ap, bq).compute_buckets()
        ({0: (-2,), 1/3: (1/3,), 1/2: (1/2, -1/2)}, {0: (1, 2)})
        """
        from collections import defaultdict

        # TODO this should probably be cached somewhere
        abuckets = defaultdict(tuple)
        bbuckets = defaultdict(tuple)

        # NOTE the new Mod object does so much canonization that we can ignore
        #      o(ab)buckets.

        for params, bucket in [(self.ap, abuckets), (self.bq, bbuckets)]:
            for p in params:
                bucket[Mod(p, 1)] += (p, )

        return dict(abuckets), dict(bbuckets)

    def build_invariants(self):
        """
        Compute the invariant vector of (``ap``, ``bq``).

        The invariant vector is:
            (gamma, ((s1, n1), ..., (sk, nk)), ((t1, m1), ..., (tr, mr)))
        where gamma is the number of integer a < 0,
              s1 < ... < sk
              nl is the number of parameters a_i congruent to sl mod 1
              t1 < ... < tr
              ml is the number of parameters b_i congruent to tl mod 1

        If the index pair contains parameters, then this is not truly an
        invariant, since the parameters cannot be sorted uniquely mod1.

        >>> from sympy.simplify.hyperexpand import IndexPair
        >>> from sympy import S
        >>> ap = (S(1)/2, S(1)/3, S(-1)/2, -2)
        >>> bq = (1, 2)

        Here gamma = 1,
             k = 3, s1 = 0, s2 = 1/3, s3 = 1/2
                    n1 = 1, n2 = 1,   n2 = 2
             r = 1, t1 = 0
                    m1 = 2:
        >>> IndexPair(ap, bq).build_invariants()
        (1, ((0, 1), (1/3, 1), (1/2, 2)), ((0, 2),))
        """
        abuckets, bbuckets = self.compute_buckets()

        gamma = 0
        if S(0) in abuckets:
            gamma = len(filter(lambda x: x < 0, abuckets[S(0)]))

        def tr(bucket):
            bucket = bucket.items()
            if not any(isinstance(x[0], Mod) for x in bucket):
                bucket.sort(key=lambda x: x[0])
            bucket = tuple(map(lambda x: (x[0], len(x[1])), bucket))
            return bucket

        return (gamma, tr(abuckets), tr(bbuckets))

    def difficulty(self, ip):
        """ Estimate how many steps it takes to reach ``ip`` from self.
            Return -1 if impossible. """
        oabuckets, obbuckets = self.compute_buckets()
        abuckets, bbuckets = ip.compute_buckets(oabuckets, obbuckets)

        gt0 = lambda x: (x > 0) is True
        if S(0) in abuckets and (not S(0) in oabuckets or
             len(filter(gt0, abuckets[S(0)])) != \
             len(filter(gt0, oabuckets[S(0)]))):
            return -1

        diff = 0
        for bucket, obucket in [(abuckets, oabuckets), (bbuckets, obbuckets)]:
            for mod in set(bucket.keys() + obucket.keys()):
                if (not mod in bucket) or (not mod in obucket) \
                   or len(bucket[mod]) != len(obucket[mod]):
                    return -1
                l1 = list(bucket[mod])
                l2 = list(obucket[mod])
                l1.sort()
                l2.sort()
                for i, j in zip(l1, l2):
                    diff += abs(i - j)

        return diff


class IndexQuadruple(object):
    """ Holds a quadruple of indices. """

    def __init__(self, an, ap, bm, bq):
        from sympy import expand, Tuple

        def tr(l):
            return Tuple(*map(expand, l))
        self.an = tr(an)
        self.ap = tr(ap)
        self.bm = tr(bm)
        self.bq = tr(bq)

    def compute_buckets(self):
        """
        Compute buckets for the fours sets of parameters.

        We guarantee that any two equal Mod objects returned are actually the
        same, and that the buckets are sorted by real part (an and bq
        descendending, bm and ap ascending).

        >>> from sympy.simplify.hyperexpand import IndexQuadruple
        >>> from sympy.abc import y
        >>> from sympy import S
        >>> a, b = [1, 3, 2, S(3)/2], [1 + y, y, 2, y + 3]
        >>> IndexQuadruple(a, b, [2], [y]).compute_buckets()
        ({0: [3, 2, 1], 1/2: [3/2]},
        {0: [2], Mod(y, 1): [y, y + 1, y + 3]}, {0: [2]}, {Mod(y, 1): [y]})

        """
        from collections import defaultdict
        dicts = pan, pap, pbm, pbq = defaultdict(list), defaultdict(list), \
                                     defaultdict(list), defaultdict(list)
        for dic, lis in zip(dicts, (self.an, self.ap, self.bm, self.bq)):
            for x in lis:
                dic[Mod(x, 1)].append(x)

        for dic, flip in zip(dicts, (True, False, False, True)):
            for m, items in dic.iteritems():
                x0 = items[0]
                items.sort(key=lambda x: x - x0, reverse=flip)
                dic[m] = items

        return tuple([dict(w) for w in dicts])

    @property
    def signature(self):
        return (len(self.an), len(self.ap), len(self.bm), len(self.bq))

    def __repr__(self):
        return 'IndexQuadruple(%s, %s, %s, %s)' % (self.an, self.ap,
                                                   self.bm, self.bq)

# Dummy variable.
_x = Dummy('x')


class Formula(object):
    """
    This class represents hypergeometric formulae.

    Its data members are:
    - z, the argument
    - closed_form, the closed form expression
    - symbols, the free symbols (parameters) in the formula
    - indices, the parameters
    - B, C, M (see _compute_basis)
    - lcms, a dictionary which maps symbol -> lcm of denominators
    - isolation, a dictonary which maps symbol -> (num, coeff) pairs

    >>> from sympy.abc import a, b, z
    >>> from sympy.simplify.hyperexpand import Formula
    >>> f = Formula((a/2, a/3 + b, (1+a)/2), (a, b, (a+b)/7), z, None, [a, b])

    The lcm of all denominators of coefficients of a is 2*3*7
    >>> f.lcms[a]
    42

    for b it is just 7:
    >>> f.lcms[b]
    7

    We can isolate a in the (1+a)/2 term, with denominator 2:
    >>> f.isolation[a]
    (2, 2, 1)

    b is isolated in the b term, with coefficient one:
    >>> f.isolation[b]
    (4, 1, 1)
    """

    def _compute_basis(self, closed_form):
        """
        Compute a set of functions B=(f1, ..., fn), a nxn matrix M
        and a 1xn matrix C such that:
           closed_form = C B
           z d/dz B = M B.
        """
        from sympy.matrices import Matrix, eye, zeros

        afactors = map(lambda a: _x + a, self.indices.ap)
        bfactors = map(lambda b: _x + b - 1, self.indices.bq)
        expr = _x*Mul(*bfactors) - self.z*Mul(*afactors)
        poly = Poly(expr, _x)

        n = poly.degree() - 1
        b = [closed_form]
        for _ in xrange(n):
            b.append(self.z*b[-1].diff(self.z))

        self.B = Matrix(b)
        self.C = Matrix([[1] + [0]*n])

        m = eye(n)
        m = m.col_insert(0, zeros(n, 1))
        l = poly.all_coeffs()[1:]
        l.reverse()
        self.M = m.row_insert(n, -Matrix([l])/poly.all_coeffs()[0])

    def __init__(self, ap, bq, z, res, symbols, B=None, C=None, M=None):
        ap = Tuple(*map(expand, ap))
        bq = Tuple(*map(expand, bq))
        z = sympify(z)
        res = sympify(res)
        symbols = filter(lambda x: ap.has(x) or bq.has(x), sympify(symbols))

        self.z = z
        self.symbols = symbols
        self.B = B
        self.C = C
        self.M = M

        params = list(ap) + list(bq)
        lcms = {}
        isolation = {}
        for a in symbols:
            from sympy import ilcm
            l = 1
            isolating = []
            others = list(symbols[:])
            others.remove(a)
            i = 0
            for p in params:
                if p.has(a):
                    c, m = None, None
                    if p.is_Add:
                        c, m = p.as_independent(a)[1].as_coeff_mul(a)
                    else:
                        c, m = p.as_coeff_mul(a)
                    if m != (a, ) or not c.is_Rational:
                        raise NotImplementedError('?')
                    l = ilcm(l, c.q)

                    if not p.has(*others):
                        isolating.append((i, c.q, c.p))
                lcms[a] = l
                i += 1
            if len(isolating) == 0:
                raise NotImplementedError('parameter is not isolated')
            isolating.sort(key=lambda x: x[1])
            isolating.sort(key=lambda x: -x[2])
            isolation[a] = isolating[-1]

        self.lcms = lcms
        self.isolation = isolation

        self.indices = IndexPair(ap, bq)

        # TODO with symbolic parameters, it could be advantageous
        #      (for prettier answers) to compute a basis only *after*
        #      instantiation
        if res is not None:
            self._compute_basis(res)

    @property
    def closed_form(self):
        return (self.C*self.B)[0]

    def find_instantiations(self, ip):
        """
        Try to find instantiations of the free symbols that match
        ``ip.ap``, ``ip.bq``. Return the instantiated formulae as a list.
        Note that the returned instantiations need not actually match,
        or be valid!
        """
        ap = ip.ap
        bq = ip.bq
        if len(ap) != len(self.indices.ap) or len(bq) != len(self.indices.bq):
            raise TypeError('Cannot instantiate other number of parameters')

        from sympy import solve
        from sympy.core.compatibility import permutations, product
        res = []
        our_params = list(self.indices.ap) + list(self.indices.bq)
        for na in permutations(ap):
            for nb in permutations(bq):
                all_params = list(na) + list(nb)
                repl = {}
                for a in self.symbols:
                    i, d, _ = self.isolation[a]
                    repl[a] = (solve(our_params[i] - all_params[i], a)[0], d)
                for change in product(*[(-1, 0, 1)]*len(self.symbols)):
                    rep = {}
                    for i, a in zip(change, repl.keys()):
                        rep[a] = repl[a][0] + i*repl[a][1]
                    res.append(Formula(self.indices.ap.subs(rep),
                                       self.indices.bq.subs(rep),
                                       self.z, None, [], self.B.subs(rep),
                                       self.C.subs(rep), self.M.subs(rep)))
                # if say a = -1/2, and there is 2*a in the formula, then
                # there will be a negative integer. But this origin is also
                # reachable from a = 1/2 ...
                # So throw this in as well.
                # The code is not as general as it could be, but good enough.
                if len(self.symbols) == 1:
                    a = self.symbols[0]
                    aval, d = repl[a]
                    if aval < 0 and d == 1:
                        from sympy import ceiling
                        aval -= ceiling(aval) - 1
                        res.append(Formula(self.indices.ap.subs(a, aval),
                                           self.indices.bq.subs(a, aval),
                                       self.z, None, [], self.B.subs(a, aval),
                                       self.C.subs(rep), self.M.subs(a, aval)))
        return res

    def is_suitable(self):
        """
        Decide if ``self`` is a suitable origin.

        >>> from sympy.simplify.hyperexpand import Formula
        >>> from sympy import S

        If ai - bq in Z and bq >= ai this is fine:
        >>> Formula((S(1)/2,), (S(3)/2,), None, None, []).is_suitable()
        True

        but ai = bq is not:
        >>> Formula((S(1)/2,), (S(1)/2,), None, None, []).is_suitable()
        False

        and ai > bq is not either:
        >>> Formula((S(1)/2,), (-S(1)/2,), None, None, []).is_suitable()
        False

        None of the bj can be a non-positive integer:
        >>> Formula((S(1)/2,), (0,), None, None, []).is_suitable()
        False
        >>> Formula((S(1)/2,), (-1, 1,), None, None, []).is_suitable()
        False

        None of the ai can be zero:
        >>> Formula((S(1)/2, 0), (1,), None, None, []).is_suitable()
        False


        More complicated examples:
        >>> Formula((S(1)/2, 1), (2, -S(2)/3), None, None, []).is_suitable()
        True
        >>> a, b = (S(1)/2, 1), (2, -S(2)/3, S(3)/2)
        >>> Formula(a, b, None, None, []).is_suitable()
        True
        """
        from sympy import oo, zoo
        if len(self.symbols) > 0:
            return None
        for a in self.indices.ap:
            for b in self.indices.bq:
                if (a - b).is_integer and not a < b:
                    return False
        for a in self.indices.ap:
            if a == 0:
                return False
        for b in self.indices.bq:
            if b <= 0 and b.is_integer:
                return False
        for e in [self.B, self.M, self.C]:
            if e is None:
                continue
            if e.has(S.NaN) or e.has(oo) or e.has(-oo) or e.has(zoo):
                return False
        return True


class FormulaCollection(object):
    """ A collection of formulae to use as origins. """

    def __init__(self):
        """ Doing this globally at module init time is a pain ... """
        self.symbolic_formulae = {}
        self.concrete_formulae = {}
        self.formulae = []

        add_formulae(self.formulae)

        # Now process the formulae into a helpful form.
        # These dicts are indexed by (p, q).

        for f in self.formulae:
            sizes = f.indices.sizes
            if len(f.symbols) > 0:
                self.symbolic_formulae.setdefault(sizes, []).append(f)
            else:
                inv = f.indices.build_invariants()
                self.concrete_formulae.setdefault(sizes, {})[inv] = f

    def lookup_origin(self, ip):
        """
        Given the suitable parameters ``ip.ap``, ``ip.bq``, try to find an
        origin in our knowledge base.

        >>> from sympy.simplify.hyperexpand import FormulaCollection, IndexPair
        >>> f = FormulaCollection()
        >>> f.lookup_origin(IndexPair((), ())).closed_form
        exp(_z)
        >>> f.lookup_origin(IndexPair([1], ())).closed_form
        HyperRep_power1(-1, _z)

        >>> from sympy import S
        >>> i = IndexPair([S('1/4'), S('3/4 + 4')], [S.Half])
        >>> f.lookup_origin(i).closed_form
        HyperRep_sqrts1(-17/4, _z)
        """
        inv = ip.build_invariants()
        sizes = ip.sizes
        if sizes in self.concrete_formulae and \
           inv in self.concrete_formulae[sizes]:
            return self.concrete_formulae[sizes][inv]

        # We don't have a concrete formula. Try to instantiate.
        if not sizes in self.symbolic_formulae:
            return None # Too bad...

        possible = []
        for f in self.symbolic_formulae[sizes]:
            l = f.find_instantiations(ip)
            for f2 in l:
                if not f2.is_suitable():
                    continue
                diff = f2.indices.difficulty(ip)
                if diff != -1:
                    possible.append((diff, f2))

        if not possible:
            # Give up.
            return None

        # find the nearest origin
        possible.sort(key=lambda x: x[0])
        return possible[0][1]


class MeijerFormula(object):
    """
    This class represents a Meijer G-function formula.

    Its data members are:
    - z, the argument
    - symbols, the free symbols (parameters) in the formula
    - indices, the parameters
    - B, C, M (c/f ordinary Formula)
    """

    def __init__(self, an, ap, bm, bq, z, symbols, B, C, M, matcher):
        an, ap, bm, bq = [Tuple(*map(expand, w)) for w in [an, ap, bm, bq]]
        self.indices = IndexQuadruple(an, ap, bm, bq)
        self.z = z
        self.symbols = symbols
        self._matcher = matcher
        self.B = B
        self.C = C
        self.M = M

    @property
    def closed_form(self):
        return (self.C*self.B)[0]

    def try_instantiate(self, iq):
        """
        Try to instantiate the current formula to (almost) match iq.
        This uses the _matcher passed on init.
        """
        if iq.signature != self.indices.signature:
            return None
        res = self._matcher(iq)
        if res is not None:
            subs, niq = res
            return MeijerFormula(niq.an, niq.ap, niq.bm, niq.bq,
                                 self.z, [],
                                 self.B.subs(subs), self.C.subs(subs),
                                 self.M.subs(subs), None)


class MeijerFormulaCollection(object):
    """
    This class holds a collection of meijer g formulae.
    """

    def __init__(self):
        from collections import defaultdict
        formulae = []
        add_meijerg_formulae(formulae)
        self.formulae = defaultdict(list)
        for formula in formulae:
            self.formulae[formula.indices.signature].append(formula)
        self.formulae = dict(self.formulae)

    def lookup_origin(self, iq):
        """ Try to find a formula that matches iq. """
        if not iq.signature in self.formulae:
            return None
        for formula in self.formulae[iq.signature]:
            res = formula.try_instantiate(iq)
            if res is not None:
                return res


class Operator(object):
    """
    Base class for operators to be applied to our functions.

    These operators are differential operators. They are by convention
    expressed in the variable D = z*d/dz (although this base class does
    not actually care).
    Note that when the operator is applied to an object, we typically do
    *not* blindly differentiate but instead use a different representation
    of the z*d/dz operator (see make_derivative_operator).

    To subclass from this, define a __init__ method that initalises a
    self._poly variable. This variable stores a polynomial. By convention
    the generator is z*d/dz, and acts to the right of all coefficients.

    Thus this poly
        x**2 + 2*z*x + 1
    represents the differential operator
        (z*d/dz)**2 + 2*z**2*d/dz.

    This class is used only in the implementation of the hypergeometric
    function expansion algorithm.
    """

    def apply(self, obj, op):
        """
        Apply ``self`` to the object ``obj``, where the generator is ``op``.

        >>> from sympy.simplify.hyperexpand import Operator
        >>> from sympy.polys.polytools import Poly
        >>> from sympy.abc import x, y, z
        >>> op = Operator()
        >>> op._poly = Poly(x**2 + z*x + y, x)
        >>> op.apply(z**7, lambda f: f.diff(z))
        y*z**7 + 7*z**7 + 42*z**5
        """
        coeffs = self._poly.all_coeffs()
        coeffs.reverse()
        diffs = [obj]
        for c in coeffs[1:]:
            diffs.append(op(diffs[-1]))
        r = coeffs[0]*diffs[0]
        for c, d in zip(coeffs[1:], diffs[1:]):
            r += c*d
        return r


class MultOperator(Operator):
    """ Simply multiply by a "constant" """

    def __init__(self, p):
        self._poly = Poly(p, _x)


class ShiftA(Operator):
    """ Increment an upper index. """

    def __init__(self, ai):
        ai = sympify(ai)
        if ai == 0:
            raise ValueError('Cannot increment zero upper index.')
        self._poly = Poly(_x/ai + 1, _x)

    def __str__(self):
        return '<Increment upper %s.>' % (1/self._poly.all_coeffs()[0])


class ShiftB(Operator):
    """ Decrement a lower index. """

    def __init__(self, bi):
        bi = sympify(bi)
        if bi == 1:
            raise ValueError('Cannot decrement unit lower index.')
        self._poly = Poly(_x/(bi - 1) + 1, _x)

    def __str__(self):
        return '<Decrement lower %s.>' % (1/self._poly.all_coeffs()[0] + 1)


class UnShiftA(Operator):
    """ Decrement an upper index. """

    def __init__(self, ap, bq, i, z):
        """ Note: i counts from zero! """
        ap, bq, i = map(sympify, [ap, bq, i])

        self._ap = ap
        self._bq = bq
        self._i = i

        ap = list(ap)
        bq = list(bq)
        ai = ap.pop(i) - 1

        if ai == 0:
            raise ValueError('Cannot decrement unit upper index.')

        m = Poly(z*ai, _x)
        for a in ap:
            m *= Poly(_x + a, _x)
        #print m

        A = Dummy('A')
        n = D = Poly(ai*A - ai, A)
        for b in bq:
            n *= (D + b - 1)
        #print n

        b0 = -n.nth(0)
        if b0 == 0:
            raise ValueError('Cannot decrement upper index: ' \
                               'cancels with lower')
        #print b0

        n = Poly(Poly(n.all_coeffs()[:-1], A).as_expr().subs(A, _x/ai + 1), _x)

        self._poly = Poly((n - m)/b0, _x)

    def __str__(self):
        return '<Decrement upper index #%s of %s, %s.>' % (self._i,
                                                        self._ap, self._bq)


class UnShiftB(Operator):
    """ Increment a lower index. """

    def __init__(self, ap, bq, i, z):
        """ Note: i counts from zero! """
        ap, bq, i = map(sympify, [ap, bq, i])

        self._ap = ap
        self._bq = bq
        self._i = i

        ap = list(ap)
        bq = list(bq)
        bi = bq.pop(i) + 1

        if bi == 0:
            raise ValueError('Cannot increment -1 lower index.')

        m = Poly(_x*(bi-1), _x)
        for b in bq:
            m *= Poly(_x + b - 1, _x)
        #print m

        B = Dummy('B')
        D = Poly((bi-1)*B - bi + 1, B)
        n = Poly(z, B)
        for a in ap:
            n *= (D + a)
        #print n

        b0 = n.nth(0)
        #print b0
        if b0 == 0:
            raise ValueError('Cannot increment index: ' \
                               'cancels with upper')
        #print b0

        n = Poly(Poly(n.all_coeffs()[:-1], B).as_expr().subs(
            B, _x/(bi-1) + 1), _x)
        #print n

        self._poly = Poly((m-n)/b0, _x)

    def __str__(self):
        return '<Increment lower index #%s of %s, %s.>' % (self._i,
                                                        self._ap, self._bq)


class MeijerShiftA(Operator):
    """ Increment an upper b index. """

    def __init__(self, bi):
        bi = sympify(bi)
        self._poly = Poly(bi - _x, _x)

    def __str__(self):
        return '<Increment upper b=%s.>' % (self._poly.all_coeffs()[1])


class MeijerShiftB(Operator):
    """ Decrement an upper a index. """

    def __init__(self, bi):
        bi = sympify(bi)
        self._poly = Poly(1 - bi + _x, _x)

    def __str__(self):
        return '<Decrement upper a=%s.>' % (1 - self._poly.all_coeffs()[1])


class MeijerShiftC(Operator):
    """ Increment a lower b index. """

    def __init__(self, bi):
        bi = sympify(bi)
        self._poly = Poly(-bi + _x, _x)

    def __str__(self):
        return '<Increment lower b=%s.>' % (-self._poly.all_coeffs()[1])


class MeijerShiftD(Operator):
    """ Decrement a lower a index. """

    def __init__(self, bi):
        bi = sympify(bi)
        self._poly = Poly(bi - 1 - _x, _x)

    def __str__(self):
        return '<Decrement lower a=%s.>' % (self._poly.all_coeffs()[1] + 1)


class MeijerUnShiftA(Operator):
    """ Decrement an upper b index. """

    def __init__(self, an, ap, bm, bq, i, z):
        """ Note: i counts from zero! """
        an, ap, bm, bq, i = map(sympify, [an, ap, bm, bq, i])

        self._an = an
        self._ap = ap
        self._bm = bm
        self._bq = bq
        self._i = i

        an = list(an)
        ap = list(ap)
        bm = list(bm)
        bq = list(bq)
        bi = bm.pop(i) - 1

        m = Poly(1, _x)
        for b in bm:
            m *= Poly(b - _x, _x)
        for b in bq:
            m *= Poly(_x - b, _x)
        #print m

        A = Dummy('A')
        D = Poly(bi - A, A)
        n = Poly(z, A)
        for a in an:
            n *= (D + 1 - a)
        for a in ap:
            n *= (-D + a - 1)
        #print n

        b0 = n.nth(0)
        #print b0
        if b0 == 0:
            raise ValueError('Cannot decrement upper b index (cancels)')
        #print b0

        n = Poly(Poly(n.all_coeffs()[:-1], A).as_expr().subs(A, bi - _x), _x)
        #print n

        self._poly = Poly((m-n)/b0, _x)

    def __str__(self):
        return '<Decrement upper b index #%s of %s, %s, %s, %s.>' % (self._i,
                                      self._an, self._ap, self._bm, self._bq)


class MeijerUnShiftB(Operator):
    """ Increment an upper a index. """

    def __init__(self, an, ap, bm, bq, i, z):
        """ Note: i counts from zero! """
        an, ap, bm, bq, i = map(sympify, [an, ap, bm, bq, i])

        self._an = an
        self._ap = ap
        self._bm = bm
        self._bq = bq
        self._i = i

        an = list(an)
        ap = list(ap)
        bm = list(bm)
        bq = list(bq)
        ai = an.pop(i) + 1

        m = Poly(z, _x)
        for a in an:
            m *= Poly(1 - a + _x, _x)
        for a in ap:
            m *= Poly(a - 1 - _x, _x)
        #print m

        B = Dummy('B')
        D = Poly(B + ai - 1, B)
        n = Poly(1, B)
        for b in bm:
            n *= (-D + b)
        for b in bq:
            n *= (D - b)
        #print n

        b0 = n.nth(0)
        #print b0
        if b0 == 0:
            raise ValueError('Cannot increment upper a index (cancels)')
        #print b0

        n = Poly(Poly(n.all_coeffs()[:-1], B).as_expr().subs(
            B, 1 - ai + _x), _x)
        #print n

        self._poly = Poly((m-n)/b0, _x)

    def __str__(self):
        return '<Increment upper a index #%s of %s, %s, %s, %s.>' % (self._i,
                                      self._an, self._ap, self._bm, self._bq)


class MeijerUnShiftC(Operator):
    """ Decrement a lower b index. """
    # XXX this is "essentially" the same as MeijerUnShiftA. This "essentially"
    #     can be made rigorous using the functional equation G(1/z) = G'(z),
    #     where G' denotes a G function of slightly altered parameters.
    #     However, sorting out the details seems harder than just coding it
    #     again.

    def __init__(self, an, ap, bm, bq, i, z):
        """ Note: i counts from zero! """
        an, ap, bm, bq, i = map(sympify, [an, ap, bm, bq, i])

        self._an = an
        self._ap = ap
        self._bm = bm
        self._bq = bq
        self._i = i

        an = list(an)
        ap = list(ap)
        bm = list(bm)
        bq = list(bq)
        bi = bq.pop(i) - 1

        m = Poly(1, _x)
        for b in bm:
            m *= Poly(b - _x, _x)
        for b in bq:
            m *= Poly(_x - b, _x)
        #print m

        C = Dummy('C')
        D = Poly(bi + C, C)
        n = Poly(z, C)
        for a in an:
            n *= (D + 1 - a)
        for a in ap:
            n *= (-D + a - 1)
        #print n

        b0 = n.nth(0)
        #print b0
        if b0 == 0:
            raise ValueError('Cannot decrement lower b index (cancels)')
        #print b0

        n = Poly(Poly(n.all_coeffs()[:-1], C).as_expr().subs(C, _x - bi), _x)
        #print n

        self._poly = Poly((m-n)/b0, _x)

    def __str__(self):
        return '<Decrement lower b index #%s of %s, %s, %s, %s.>' % (self._i,
                                      self._an, self._ap, self._bm, self._bq)


class MeijerUnShiftD(Operator):
    """ Increment a lower a index. """
    # XXX This is essentially the same as MeijerUnShiftA.
    #     See comment at MeijerUnShiftC.

    def __init__(self, an, ap, bm, bq, i, z):
        """ Note: i counts from zero! """
        an, ap, bm, bq, i = map(sympify, [an, ap, bm, bq, i])

        self._an = an
        self._ap = ap
        self._bm = bm
        self._bq = bq
        self._i = i

        an = list(an)
        ap = list(ap)
        bm = list(bm)
        bq = list(bq)
        ai = ap.pop(i) + 1

        m = Poly(z, _x)
        for a in an:
            m *= Poly(1 - a + _x, _x)
        for a in ap:
            m *= Poly(a - 1 - _x, _x)
        #print m

        B = Dummy('B') # - this is the shift operator `D_I`
        D = Poly(ai - 1 - B, B)
        n = Poly(1, B)
        for b in bm:
            n *= (-D + b)
        for b in bq:
            n *= (D - b)
        #print n

        b0 = n.nth(0)
        #print b0
        if b0 == 0:
            raise ValueError('Cannot increment lower a index (cancels)')
        #print b0

        n = Poly(Poly(n.all_coeffs()[:-1], B).as_expr().subs(
            B, ai - 1 - _x), _x)
        #print n

        self._poly = Poly((m-n)/b0, _x)

    def __str__(self):
        return '<Increment lower a index #%s of %s, %s, %s, %s.>' % (self._i,
                                      self._an, self._ap, self._bm, self._bq)


class ReduceOrder(Operator):
    """ Reduce Order by cancelling an upper and a lower index. """

    def __new__(cls, ai, bj):
        """ For convenience if reduction is not possible, return None. """
        ai = sympify(ai)
        bj = sympify(bj)
        n = ai - bj
        if n < 0 or not n.is_Integer:
            return None
        if bj.is_integer and bj <= 0 and bj + n - 1 >= 0:
            return None

        self = Operator.__new__(cls)

        p = S(1)
        for k in xrange(n):
            p *= (_x + bj + k)/(bj + k)

        self._poly = Poly(p, _x)
        self._a = ai
        self._b = bj

        return self

    @classmethod
    def _meijer(cls, b, a, sign):
        """ Cancel b + sign*s and a + sign*s
            This is for meijer G functions. """
        from sympy import Add
        b = sympify(b)
        a = sympify(a)
        n = b - a
        if n < 0 or not n.is_Integer:
            return None

        self = Operator.__new__(cls)

        p = S(1)
        for k in xrange(n):
            p *= (sign*_x + a + k)

        self._poly = Poly(p, _x)
        if sign == -1:
            self._a = b
            self._b = a
        else:
            self._b = Add(1, a - 1, evaluate=False)
            self._a = Add(1, b - 1, evaluate=False)

        return self

    @classmethod
    def meijer_minus(cls, b, a):
        return cls._meijer(b, a, -1)

    @classmethod
    def meijer_plus(cls, a, b):
        return cls._meijer(1 - a, 1 - b, 1)

    def __str__(self):
        return '<Reduce order by cancelling upper %s with lower %s.>' % \
                  (self._a, self._b)


def _reduce_order(ap, bq, gen, key):
    """ Order reduction algorithm used in Hypergeometric and Meijer G """
    ap = list(ap)
    bq = list(bq)

    ap.sort(key=key)
    bq.sort(key=key)

    nap = []
    # we will edit bq in place
    operators = []
    for a in ap:
        op = None
        for i in xrange(len(bq)):
            op = gen(a, bq[i])
            if op is not None:
                bq.pop(i)
                break
        if op is None:
            nap.append(a)
        else:
            operators.append(op)

    return nap, bq, operators


def reduce_order(ip):
    """
    Given the hypergeometric parameters ``ip.ap``, ``ip.bq``,
    find a sequence of operators to reduces order as much as possible.

    Return (nip, [operators]), where applying the operators to the
    hypergeometric function specified by nip.ap, nip.bq yields ap, bq.

    Examples
    ========

    >>> from sympy.simplify.hyperexpand import reduce_order, IndexPair
    >>> reduce_order(IndexPair((1, 2), (3, 4)))
    (IndexPair((1, 2), (3, 4)), [])
    >>> reduce_order(IndexPair((1,), (1,)))
    (IndexPair((), ()), [<Reduce order by cancelling upper 1 with lower 1.>])
    >>> reduce_order(IndexPair((2, 4), (3, 3)))
    (IndexPair((2,), (3,)), [<Reduce order by cancelling
    upper 4 with lower 3.>])
    """
    nap, nbq, operators = _reduce_order(ip.ap, ip.bq, ReduceOrder, lambda x: x)

    return IndexPair(Tuple(*nap), Tuple(*nbq)), operators


def reduce_order_meijer(iq):
    """
    Given the Meijer G function parameters, ``iq.am``, ``iq.ap``, ``iq.bm``,
    ``iq.bq``, find a sequence of operators that reduces order as much as
    possible.

    Return niq, [operators].

    Examples
    ========

    >>> from sympy.simplify.hyperexpand import (reduce_order_meijer,
    ...                                         IndexQuadruple)
    >>> reduce_order_meijer(IndexQuadruple([3, 4], [5, 6], [3, 4], [1, 2]))[0]
    IndexQuadruple((4, 3), (5, 6), (3, 4), (2, 1))
    >>> reduce_order_meijer(IndexQuadruple([3, 4], [5, 6], [3, 4], [1, 8]))[0]
    IndexQuadruple((3,), (5, 6), (3, 4), (1,))
    >>> reduce_order_meijer(IndexQuadruple([3, 4], [5, 6], [7, 5], [1, 5]))[0]
    IndexQuadruple((3,), (), (), (1,))
    >>> reduce_order_meijer(IndexQuadruple([3, 4], [5, 6], [7, 5], [5, 3]))[0]
    IndexQuadruple((), (), (), ())
    """

    nan, nbq, ops1 = _reduce_order(iq.an, iq.bq, ReduceOrder.meijer_plus,
                                   lambda x: -x)
    nbm, nap, ops2 = _reduce_order(iq.bm, iq.ap, ReduceOrder.meijer_minus,
                                   lambda x: x)

    return IndexQuadruple(*[Tuple(*w) for w in [nan, nap, nbm, nbq]]), \
           ops1 + ops2


def make_derivative_operator(M, z):
    """ Create a derivative operator, to be passed to Operator.apply. """
    from sympy import poly

    def doit(C):
        r = z*C.diff(z) + C*M
        r = r.applyfunc(make_simp(z))
        return r
    return doit


def apply_operators(obj, ops, op):
    """
    Apply the list of operators ``ops`` to object ``obj``, substituting
    ``op`` for the generator.
    """
    res = obj
    for o in reversed(ops):
        res = o.apply(res, op)
    return res


def devise_plan(ip, nip, z):
    """
    Devise a plan (consisting of shift and un-shift operators) to be applied
    to the hypergeometric function (``nip.ap``, ``nip.bq``) to yield
    (``ip.ap``, ``ip.bq``).
    Returns a list of operators.

    >>> from sympy.simplify.hyperexpand import devise_plan, IndexPair
    >>> from sympy.abc import z

    Nothing to do:

    >>> devise_plan(IndexPair((1, 2), ()), IndexPair((1, 2), ()), z)
    []
    >>> devise_plan(IndexPair((), (1, 2)), IndexPair((), (1, 2)), z)
    []

    Very simple plans:

    >>> devise_plan(IndexPair((2,), ()), IndexPair((1,), ()), z)
    [<Increment upper 1.>]
    >>> devise_plan(IndexPair((), (2,)), IndexPair((), (1,)), z)
    [<Increment lower index #0 of [], [1].>]

    Several buckets:

    >>> from sympy import S
    >>> devise_plan(IndexPair((1, S.Half), ()),
    ...             IndexPair((2, S('3/2')), ()), z) #doctest: +NORMALIZE_WHITESPACE
    [<Decrement upper index #0 of [3/2, 1], [].>,
    <Decrement upper index #0 of [2, 3/2], [].>]

    A slightly more complicated plan:

    >>> devise_plan(IndexPair((1, 3), ()), IndexPair((2, 2), ()), z)
    [<Increment upper 2.>, <Decrement upper index #0 of [2, 2], [].>]

    Another more complicated plan: (note that the ap have to be shifted first!)

    >>> devise_plan(IndexPair((1, -1), (2,)), IndexPair((3, -2), (4,)), z)
    [<Decrement lower 3.>, <Decrement lower 4.>,
    <Decrement upper index #1 of [-1, 2], [4].>,
    <Decrement upper index #1 of [-1, 3], [4].>, <Increment upper -2.>]
    """
    abuckets, bbuckets = ip.compute_buckets()
    nabuckets, nbbuckets = nip.compute_buckets(abuckets, bbuckets)

    if len(abuckets.keys()) != len(nabuckets.keys()) or \
       len(bbuckets.keys()) != len(nbbuckets.keys()):
        raise ValueError('%s not reachable from %s' % (ip, nip))

    ops = []

    def do_shifts(fro, to, inc, dec):
        ops = []
        for i in xrange(len(fro)):
            if to[i] - fro[i] > 0:
                sh = inc
                ch = 1
            else:
                sh = dec
                ch = -1

            while to[i] != fro[i]:
                ops += [sh(fro, i)]
                fro[i] += ch

        return ops

    def do_shifts_a(nal, nbk, al, aother, bother):
        """ Shift us from (nal, nbk) to (al, nbk). """
        return do_shifts(nal, al, lambda p, i: ShiftA(p[i]),
                         lambda p, i: UnShiftA(p + aother, nbk + bother, i, z))

    def do_shifts_b(nal, nbk, bk, aother, bother):
        """ Shift us from (nal, nbk) to (nal, bk). """
        return do_shifts(nbk, bk,
                         lambda p, i: UnShiftB(nal + aother, p + bother, i, z),
                         lambda p, i: ShiftB(p[i]))

    for r in sorted(abuckets.keys() + bbuckets.keys(), key=default_sort_key):
        al = ()
        nal = ()
        bk = ()
        nbk = ()
        if r in abuckets:
            al = abuckets[r]
            nal = nabuckets[r]
        if r in bbuckets:
            bk = bbuckets[r]
            nbk = nbbuckets[r]
        if len(al) != len(nal) or len(bk) != len(nbk):
            raise ValueError('%s not reachable from %s' % (
                             (ip.ap, ip.bq), (nip.ap, nip.bq)))

        al, nal, bk, nbk = [sorted(list(w), key=default_sort_key)
            for w in [al, nal, bk, nbk]]

        def others(dic, key):
            l = []
            for k, value in dic.iteritems():
                if k != key:
                    l += list(dic[k])
            return l
        aother = others(nabuckets, r)
        bother = others(nbbuckets, r)

        if len(al) == 0:
            # there can be no complications, just shift the bs as we please
            ops += do_shifts_b([], nbk, bk, aother, bother)
        elif len(bk) == 0:
            # there can be no complications, just shift the as as we please
            ops += do_shifts_a(nal, [], al, aother, bother)
        else:
            namax = nal[-1]
            amax = al[-1]

            if nbk[0] <= namax or bk[0] <= amax:
                raise ValueError('Non-suitable parameters.')

            if namax > amax:
                # we are going to shift down - first do the as, then the bs
                ops += do_shifts_a(nal, nbk, al, aother, bother)
                ops += do_shifts_b(al, nbk, bk, aother, bother)
            else:
                # we are going to shift up - first do the bs, then the as
                ops += do_shifts_b(nal, nbk, bk, aother, bother)
                ops += do_shifts_a(nal, bk, al, aother, bother)

        nabuckets[r] = al
        nbbuckets[r] = bk

    ops.reverse()
    return ops


def try_shifted_sum(ip, z):
    """ Try to recognise a hypergeometric sum that starts from k > 0. """
    from sympy.functions import rf, factorial
    abuckets, bbuckets = ip.compute_buckets()
    if not S(0) in abuckets or len(abuckets[S(0)]) != 1:
        return None
    r = abuckets[S(0)][0]
    if r <= 0:
        return None
    if not S(0) in bbuckets:
        return None
    l = list(bbuckets[S(0)])
    l.sort()
    k = l[0]
    if k <= 0:
        return None

    nap = list(ip.ap)
    nap.remove(r)
    nbq = list(ip.bq)
    nbq.remove(k)
    k -= 1
    nap = map(lambda x: x - k, nap)
    nbq = map(lambda x: x - k, nbq)

    ops = []
    for n in xrange(r - 1):
        ops.append(ShiftA(n + 1))
    ops.reverse()

    fac = factorial(k)/z**k
    for a in nap:
        fac /= rf(a, k)
    for b in nbq:
        fac *= rf(b, k)

    ops += [MultOperator(fac)]

    p = 0
    for n in xrange(k):
        m = z**n/factorial(n)
        for a in nap:
            m *= rf(a, n)
        for b in nbq:
            m /= rf(b, n)
        p += m

    return IndexPair(nap, nbq), ops, -p


def try_polynomial(ip, z):
    """ Recognise polynomial cases. Returns None if not such a case.
        Requires order to be fully reduced. """
    from sympy import oo, factorial, rf, Expr
    abuckets, bbuckets = ip.compute_buckets()
    a0 = list(abuckets.get(S(0), []))
    b0 = list(bbuckets.get(S(0), []))
    a0.sort()
    b0.sort()
    al0 = filter(lambda x: x <= 0, a0)
    bl0 = filter(lambda x: x <= 0, b0)

    if bl0:
        return oo
    if not al0:
        return None

    a = al0[-1]
    fac = 1
    res = S(1)
    for n in Tuple(*range(-a)):
        fac *= z
        fac /= n + 1
        for a in ip.ap:
            fac *= a + n
        for b in ip.bq:
            fac /= b + n
        res += fac
    return res


def try_lerchphi(nip):
    """
    Try to find an expression for IndexPair ``nip`` in terms of Lerch
    Transcendents.

    Return None if no such expression can be found.
    """
    # This is actually quite simple, and is described in Roach's paper,
    # section 18.
    # We don't need to implement the reduction to polylog here, this
    # is handled by expand_func.
    from sympy import (expand_func, lerchphi, apart, Dummy, rf, Poly, Matrix,
                       zeros, Add)

    # First we need to figure out if the summation coefficient is a rational
    # function of the summation index, and construct that rational function.
    abuckets, bbuckets = nip.compute_buckets()
    # Update all the keys in bbuckets to be the same Mod objects as abuckets.
    akeys = abuckets.keys()
    bkeys = []
    for key in bbuckets.keys():
        new = key
        for a in akeys:
            if a == key:
                new = a
                break
        bkeys += [(new, key)]
    bb = {}
    for new, key in bkeys:
        bb[new] = bbuckets[key]
    bbuckets = bb

    paired = {}
    for key, value in abuckets.items():
        if key != 0 and not key in bbuckets:
            return None
        bvalue = bbuckets.get(key, [])
        paired[key] = (list(value), list(bvalue))
        bbuckets.pop(key, None)
    if bbuckets != {}:
        return None
    if not S(0) in abuckets:
        return None
    aints, bints = paired[S(0)]
    # Account for the additional n! in denominator
    paired[S(0)] = (aints, bints + [1])

    t = Dummy('t')
    numer = S(1)
    denom = S(1)
    for key, (avalue, bvalue) in paired.items():
        if len(avalue) != len(bvalue):
            return None
        # Note that since order has been reduced fully, all the b are
        # bigger than all the a they differ from by an integer. In particular
        # if there are any negative b left, this function is not well-defined.
        for a, b in zip(avalue, bvalue):
            if a > b:
                k = a - b
                numer *= rf(b + t, k)
                denom *= rf(b, k)
            else:
                k = b - a
                numer *= rf(a, k)
                denom *= rf(a + t, k)

    # Now do a partial fraction decomposition.
    # We assemble two structures: a list monomials of pairs (a, b) representing
    # a*t**b (b a non-negative integer), and a dict terms, where
    # terms[a] = [(b, c)] means that there is a term b/(t-a)**c.
    part = apart(numer/denom, t)
    args = Add.make_args(part)
    monomials = []
    terms = {}
    for arg in args:
        numer, denom = arg.as_numer_denom()
        if not denom.has(t):
            p = Poly(numer, t)
            assert p.is_monomial
            ((b, ), a) = p.LT()
            monomials += [(a/denom, b)]
            continue
        if numer.has(t):
            raise NotImplementedError('Need partial fraction decomposition' \
                                      ' with linear denominators')
        indep, [dep] = denom.as_coeff_mul(t)
        n = 1
        if dep.is_Pow:
            n = dep.exp
            dep = dep.base
        if dep == t:
            a == 0
        elif dep.is_Add:
            a, tmp = dep.as_independent(t)
            b = 1
            if tmp != t:
                b, _ = tmp.as_independent(t)
            if dep != b*t + a:
                raise NotImplementedError('unrecognised form %s' % dep)
            a /= b
            indep *= b**n
        else:
            raise NotImplementedError('unrecognised form of partial fraction')
        terms.setdefault(a, []).append((numer/indep, n))

    # Now that we have this information, assemble our formula. All the
    # monomials yield rational functions and go into one basis element.
    # The terms[a] are related by differentiation. If the largest exponent is
    # n, we need lerchphi(z, k, a) for k = 1, 2, ..., n.
    # deriv maps a basis to its derivative, expressed as a C(z)-linear
    # combination of other basis elements.
    deriv = {}
    coeffs = {}
    z = Dummy('z')
    monomials.sort(key=lambda x: x[1])
    mon = {0: 1/(1-z)}
    if monomials:
        for k in range(monomials[-1][1]):
            mon[k + 1] = z*mon[k].diff(z)
    for a, n in monomials:
        coeffs.setdefault(S(1), []).append(a*mon[n])
    for a, l in terms.items():
        for c, k in l:
            coeffs.setdefault(lerchphi(z, k, a), []).append(c)
        l.sort(key=lambda x: x[1])
        for k in range(2, l[-1][1] + 1):
            deriv[lerchphi(z, k, a)] = [(-a, lerchphi(z, k, a)),
                                        (1, lerchphi(z, k-1, a))]
        deriv[lerchphi(z, 1, a)] = [(-a, lerchphi(z, 1, a)),
                                    (1/(1 - z), S(1))]
    trans = {}
    for n, b in enumerate([S(1)] + deriv.keys()):
        trans[b] = n
    basis = [expand_func(b) for (b, _) in sorted(trans.items(),
                                                 key=lambda x:x[1])]
    B = Matrix(basis)
    C = Matrix([[0]*len(B)])
    for b, c in coeffs.items():
        C[trans[b]] = Add(*c)
    M = zeros(len(B))
    for b, l in deriv.items():
        for c, b2 in l:
            M[trans[b], trans[b2]] = c
    return Formula(nip.ap, nip.bq, z, None, [], B, C, M)


def build_hypergeometric_formula(nip):
    """ Create a formula object representing the hypergeometric function
        Corresponding to nip. """
    # We know that no `ap` are negative integers, otherwise "detect poly"
    # would have kicked in. However, `ap` could be empty. In this case we can
    # use a different basis.
    # I'm not aware of a basis that works in all cases.
    from sympy import zeros, Dummy, Matrix, hyper, eye, Mul
    z = Dummy('z')
    if nip.ap:
        afactors = map(lambda a: _x + a, nip.ap)
        bfactors = map(lambda b: _x + b - 1, nip.bq)
        expr = _x*Mul(*bfactors) - z*Mul(*afactors)
        poly = Poly(expr, _x)
        n = poly.degree()
        basis = []
        M = zeros(n)
        for k in xrange(n):
            a = nip.ap[0] + k
            basis += [hyper([a] + list(nip.ap[1:]), nip.bq, z)]
            if k < n - 1:
                M[k, k] = -a
                M[k, k + 1] = a
        B = Matrix(basis)
        C = Matrix([[1] + [0]*(n - 1)])
        derivs = [eye(n)]
        for k in xrange(n):
            derivs.append(M*derivs[k])
        l = poly.all_coeffs()
        l.reverse()
        res = [0]*n
        for k, c in enumerate(l):
            for r, d in enumerate(C*derivs[k]):
                res[r] += c*d
        for k, c in enumerate(res):
            M[n-1, k] = -c/derivs[n-1][0, n-1]/poly.all_coeffs()[0]
        return Formula(nip.ap, nip.bq, z, None, [], B, C, M)
    else:
        # Since there are no `ap`, none of the `bq` can be non-positive
        # integers.
        basis = []
        bq = list(nip.bq[:])
        for i in range(len(bq)):
            basis += [hyper([], bq, z)]
            bq[i] += 1
        basis += [hyper([], bq, z)]
        B = Matrix(basis)
        n = len(B)
        C = Matrix([[1] + [0]*(n - 1)])
        M = zeros(n)
        M[0, n - 1] = z/Mul(*nip.bq)
        for k in range(1, n):
            M[k, k - 1] = nip.bq[k - 1]
            M[k, k] = -nip.bq[k - 1]
        return Formula(nip.ap, nip.bq, z, None, [], B, C, M)


def hyperexpand_special(ap, bq, z):
    """
    Try to find a closed-form expression for hyper(ap, bq, z), where ``z``
    is supposed to be a "special" value, e.g. 1.

    This function tries various of the classical summation formulae
    (Gauss, Saalschuetz, etc).
    """
    # This code is very ad-hoc. There are many clever algorithms
    # (notably Zeilberger's) related to this problem.
    # For now we just want a few simple cases to work.
    from sympy import gamma, simplify, cos, unpolarify, pi
    p, q = len(ap), len(bq)
    z_ = z
    z = unpolarify(z)
    if z == 0:
        return S.Zero
    if p == 2 and q == 1:
        # 2F1
        a, b, c = ap + bq
        if z == 1:
            # Gauss
            return gamma(c - a - b)*gamma(c)/gamma(c - a)/gamma(c - b)
        if z == -1 and simplify(b - a + c) == 1:
            b, a = a, b
        if z == -1 and simplify(a - b + c) == 1:
            # Kummer
            if b.is_integer and b < 0:
                return 2*cos(pi*b/2)*gamma(-b)*gamma(b - a + 1) \
                       /gamma(-b/2)/gamma(b/2 - a + 1)
            else:
                return gamma(b/2 + 1)*gamma(b - a + 1) \
                       /gamma(b + 1)/gamma(b/2 - a + 1)
    # TODO tons of more formulae
    #      investigate what algorithms exist
    return hyper(ap, bq, z_)

_collection = None


@_timeit
def _hyperexpand(ip, z, ops0=[], z0=Dummy('z0'), premult=1, prem=0,
                 rewrite='default'):
    """
    Try to find an expression for the hypergeometric function
    ``ip.ap``, ``ip.bq``.

    The result is expressed in terms of a dummy variable z0. Then it
    is multiplied by premult. Then ops0 is applied.
    premult must be a*z**prem for some a independent of z.
    """
    from sympy.simplify import powdenest, simplify, polarify, unpolarify
    z = polarify(z, subs=False)
    if rewrite == 'default':
        rewrite = 'nonrepsmall'

    def carryout_plan(f, ops):
        C = apply_operators(f.C.subs(f.z, z0), ops,
                            make_derivative_operator(f.M.subs(f.z, z0), z0))
        from sympy import eye
        C = apply_operators(C, ops0,
                            make_derivative_operator(f.M.subs(f.z, z0)
                                         + prem*eye(f.M.shape[0]), z0))

        if premult == 1:
            C = C.applyfunc(make_simp(z0))
        r = C*f.B.subs(f.z, z0)*premult
        res = r[0].subs(z0, z)
        if rewrite:
            res = res.rewrite(rewrite)
        return res

    # TODO
    # The following would be possible:
    # *) PFD Duplication (see Kelly Roach's paper)
    # *) In a similar spirit, try_lerchphi() can be generalised considerably.

    global _collection
    if _collection is None:
        _collection = FormulaCollection()

    debug('Trying to expand hypergeometric function corresponding to', ip)

    # First reduce order as much as possible.
    nip, ops = reduce_order(ip)
    if ops:
        debug('  Reduced order to', nip)
    else:
        debug('  Could not reduce order.')

    # Now try polynomial cases
    res = try_polynomial(nip, z0)
    if res is not None:
        debug('  Recognised polynomial.')
        p = apply_operators(res, ops, lambda f: z0*f.diff(z0))
        p = apply_operators(p*premult, ops0, lambda f: z0*f.diff(z0))
        return unpolarify(simplify(p).subs(z0, z))

    # Try to recognise a shifted sum.
    p = S(0)
    res = try_shifted_sum(nip, z0)
    if res != None:
        nip, nops, p = res
        debug('  Recognised shifted sum, reducerd order to', nip)
        ops += nops

    # apply the plan for poly
    p = apply_operators(p, ops, lambda f: z0*f.diff(z0))
    p = apply_operators(p*premult, ops0, lambda f: z0*f.diff(z0))
    p = simplify(p).subs(z0, z)

    # Try special expansions early.
    if unpolarify(z) in [1, -1] and (len(nip.ap), len(nip.bq)) == (2, 1):
        f = build_hypergeometric_formula(nip)
        r = carryout_plan(f, ops).replace(hyper, hyperexpand_special)
        if not r.has(hyper):
            return r + p

    # Try to find a formula in our collection
    f = _collection.lookup_origin(nip)

    # Now try a lerch phi formula
    if f is None:
        f = try_lerchphi(nip)

    if f is None:
        debug('  Could not find an origin.',
              'Will return answer in terms of '+
              'simpler hypergeometric functions.')
        f = build_hypergeometric_formula(nip)

    debug('  Found an origin:', f.closed_form, f.indices)

    # We need to find the operators that convert f into (nap, nbq).
    ops += devise_plan(nip, f.indices, z0)

    # Now carry out the plan.
    r = carryout_plan(f, ops) + p

    return powdenest(r, polar=True).replace(hyper, hyperexpand_special)


def devise_plan_meijer(fro, to, z):
    """
    Find a sequence of operators to convert index quadruple ``fro`` into
    index quadruple ``to``. It is assumed that fro and to have the same
    signatures, and that in fact any corresponding pair of parameters differs
    by integers, and a direct path is possible. I.e. if there are parameters
       a1 b1 c1  and a2 b2 c2
    it is assumed that a1 can be shifted to a2, etc.
    The only thing this routine determines is the order of shifts to apply,
    nothing clever will be tried.
    It is also assumed that fro is suitable.

    >>> from sympy.simplify.hyperexpand import (devise_plan_meijer,
    ...                                         IndexQuadruple)
    >>> from sympy.abc import z

    Empty plan:

    >>> devise_plan_meijer(IndexQuadruple([1], [2], [3], [4]),
    ...                    IndexQuadruple([1], [2], [3], [4]), z)
    []

    Very simple plans:

    >>> devise_plan_meijer(IndexQuadruple([0], [], [], []),
    ...                    IndexQuadruple([1], [], [], []), z)
    [<Increment upper a index #0 of [0], [], [], [].>]
    >>> devise_plan_meijer(IndexQuadruple([0], [], [], []),
    ...                    IndexQuadruple([-1], [], [], []), z)
    [<Decrement upper a=0.>]
    >>> devise_plan_meijer(IndexQuadruple([], [1], [], []),
    ...                    IndexQuadruple([], [2], [], []), z)
    [<Increment lower a index #0 of [], [1], [], [].>]

    Slightly more complicated plans:

    >>> devise_plan_meijer(IndexQuadruple([0], [], [], []),
    ...                    IndexQuadruple([2], [], [], []), z)
    [<Increment upper a index #0 of [1], [], [], [].>,
    <Increment upper a index #0 of [0], [], [], [].>]
    >>> devise_plan_meijer(IndexQuadruple([0], [], [0], []),
    ...                    IndexQuadruple([-1], [], [1], []), z)
    [<Increment upper b=0.>, <Decrement upper a=0.>]

    Order matters:

    >>> devise_plan_meijer(IndexQuadruple([0], [], [0], []),
    ...                    IndexQuadruple([1], [], [1], []), z)
    [<Increment upper a index #0 of [0], [], [1], [].>, <Increment upper b=0.>]
    """
    # TODO for now, we use the following simple heuristic: inverse-shift
    #      when possible, shift otherwise. Give up if we cannot make progress.

    def try_shift(f, t, shifter, diff, counter):
        """ Try to apply ``shifter`` in order to bring some element in ``f``
            nearer to its counterpart in ``to``. ``diff`` is +/- 1 and
            determines the effect of ``shifter``. Counter is a list of elements
            blocking the shift.

            Return an operator if change was possible, else None.
        """
        for idx, (a, b) in enumerate(zip(f, t)):
            if (
                (a - b).is_integer and (b - a)/diff > 0 and
                all(a != x for x in counter)):
                sh = shifter(idx)
                f[idx] += diff
                return sh
    fan = list(fro.an)
    fap = list(fro.ap)
    fbm = list(fro.bm)
    fbq = list(fro.bq)
    ops = []
    change = True
    while change:
        change = False
        op = try_shift(fan, to.an,
                       lambda i: MeijerUnShiftB(fan, fap, fbm, fbq, i, z),
                       1, fbm + fbq)
        if op is not None:
            ops += [op]
            change = True
            continue
        op = try_shift(fap, to.ap,
                       lambda i: MeijerUnShiftD(fan, fap, fbm, fbq, i, z),
                       1, fbm + fbq)
        if op is not None:
            ops += [op]
            change = True
            continue
        op = try_shift(fbm, to.bm,
                       lambda i: MeijerUnShiftA(fan, fap, fbm, fbq, i, z),
                       -1, fan + fap)
        if op is not None:
            ops += [op]
            change = True
            continue
        op = try_shift(fbq, to.bq,
                       lambda i: MeijerUnShiftC(fan, fap, fbm, fbq, i, z),
                       -1, fan + fap)
        if op is not None:
            ops += [op]
            change = True
            continue
        op = try_shift(fan, to.an, lambda i: MeijerShiftB(fan[i]), -1, [])
        if op is not None:
            ops += [op]
            change = True
            continue
        op = try_shift(fap, to.ap, lambda i: MeijerShiftD(fap[i]), -1, [])
        if op is not None:
            ops += [op]
            change = True
            continue
        op = try_shift(fbm, to.bm, lambda i: MeijerShiftA(fbm[i]), 1, [])
        if op is not None:
            ops += [op]
            change = True
            continue
        op = try_shift(fbq, to.bq, lambda i: MeijerShiftC(fbq[i]), 1, [])
        if op is not None:
            ops += [op]
            change = True
            continue
    if fan != list(to.an) or fap != list(to.ap) or fbm != list(to.bm) or \
       fbq != list(to.bq):
        raise NotImplementedError('Could not devise plan.')
    ops.reverse()
    return ops

_meijercollection = None


@_timeit
def _meijergexpand(iq, z0, allow_hyper=False, rewrite='default'):
    """
    Try to find an expression for the Meijer G function specified
    by the IndexQuadruple ``iq``. If ``allow_hyper`` is True, then returning
    an expression in terms of hypergeometric functions is allowed.

    Currently this just does slater's theorem.
    """
    from sympy import hyper, Piecewise, meijerg, powdenest, re, polar_lift, oo
    global _meijercollection
    if _meijercollection is None:
        _meijercollection = MeijerFormulaCollection()
    if rewrite == 'default':
        rewrite = None

    iq_ = iq
    debug('Try to expand meijer G function corresponding to', iq)

    # We will play games with analytic continuation - rather use a fresh symbol
    z = Dummy('z')

    iq, ops = reduce_order_meijer(iq)
    if ops:
        debug('  Reduced order to', iq)
    else:
        debug('  Could not reduce order.')

    # Try to find a direct formula
    f = _meijercollection.lookup_origin(iq)
    if f is not None:
        debug('  Found a Meijer G formula:', f.indices)
        ops += devise_plan_meijer(f.indices, iq, z)

        # Now carry out the plan.
        C = apply_operators(f.C.subs(f.z, z), ops,
                            make_derivative_operator(f.M.subs(f.z, z), z))

        C = C.applyfunc(make_simp(z))
        r = C*f.B.subs(f.z, z)
        r = r[0].subs(z, z0)
        return powdenest(r, polar=True)

    debug("  Could not find a direct formula. Trying slater's theorem.")

    # TODO the following would be possible:
    # *) Paired Index Theorems
    # *) PFD Duplication
    #    (See Kelly Roach's paper for details on either.)
    #
    # TODO Also, we tend to create combinations of gamma functions that can be
    #      simplified.

    def can_do(pbm, pap):
        """ Test if slater applies. """
        for i in pbm:
            if len(pbm[i]) > 1:
                l = 0
                if i in pap:
                    l = len(pap[i])
                if l + 1 < len(pbm[i]):
                    return False
        return True

    def do_slater(an, bm, ap, bq, z, zfinal):
        from sympy import gamma, residue, factorial, rf, expand_func, \
                          polar_lift
        # zfinal is the value that will eventually be substituted for z.
        # We pass it to _hyperexpand to improve performance.
        iq = IndexQuadruple(an, bm, ap, bq)
        _, pbm, pap, _ = iq.compute_buckets()
        if not can_do(pbm, pap):
            return S(0), False

        cond = len(an) + len(ap) < len(bm) + len(bq)
        if len(an) + len(ap) == len(bm) + len(bq):
            cond = abs(z) < 1
        if cond is False:
            return S(0), False

        res = S(0)
        for m in pbm:
            if len(pbm[m]) == 1:
                bh = pbm[m][0]
                fac = 1
                bo = list(bm)
                bo.remove(bh)
                for bj in bo:
                    fac *= gamma(bj - bh)
                for aj in an:
                    fac *= gamma(1 + bh - aj)
                for bj in bq:
                    fac /= gamma(1 + bh - bj)
                for aj in ap:
                    fac /= gamma(aj - bh)
                nap = [1 + bh - a for a in list(an) + list(ap)]
                nbq = [1 + bh - b for b in list(bo) + list(bq)]

                k = polar_lift(S(-1)**(len(ap) - len(bm)))
                harg = k*zfinal
                # NOTE even though k "is" +-1, this has to be t/k instead of
                #      t*k ... we are using polar numbers for consistency!
                premult = (t/k)**bh
                hyp = _hyperexpand(IndexPair(nap, nbq), harg, ops,
                                   t, premult, bh, rewrite=None)
                res += fac * hyp
            else:
                b_ = pbm[m][0]
                ki = [bi - b_ for bi in pbm[m][1:]]
                u = len(ki)
                li = [ai - b_ for ai in pap[m][:u+1]]
                bo = list(bm)
                for b in pbm[m]:
                    bo.remove(b)
                ao = list(ap)
                for a in pap[m][:u]:
                    ao.remove(a)
                lu = li[-1]
                di = [l - k for (l, k) in zip(li, ki)]

                # We first work out the integrand:
                s = Dummy('s')
                integrand = z**s
                for b in bm:
                    integrand *= gamma(b - s)
                for a in an:
                    integrand *= gamma(1 - a + s)
                for b in bq:
                    integrand /= gamma(1 - b + s)
                for a in ap:
                    integrand /= gamma(a - s)

                # Now sum the finitely many residues:
                # XXX This speeds up some cases - is it a good idea?
                integrand = expand_func(integrand)
                for r in range(lu):
                    resid = residue(integrand, s, b_ + r)
                    resid = apply_operators(resid, ops, lambda f: z*f.diff(z))
                    res -= resid

                # Now the hypergeometric term.
                au = b_ + lu
                k = polar_lift(S(-1)**(len(ao) + len(bo) + 1))
                harg = k*zfinal
                premult = (t/k)**au
                nap = [1 + au - a for a in list(an) + list(ap)] + [1]
                nbq = [1 + au - b for b in list(bm) + list(bq)]

                hyp = _hyperexpand(IndexPair(nap, nbq), harg, ops,
                                   t, premult, au, rewrite=None)

                C = S(-1)**(lu)/factorial(lu)
                for i in range(u):
                    C *= S(-1)**di[i]/rf(lu - li[i] + 1, di[i])
                for a in an:
                    C *= gamma(1 - a + au)
                for b in bo:
                    C *= gamma(b - au)
                for a in ao:
                    C /= gamma(a - au)
                for b in bq:
                    C /= gamma(1 - b + au)

                res += C*hyp

        return res, cond

    t = Dummy('t')
    slater1, cond1 = do_slater(iq.an, iq.bm, iq.ap, iq.bq, z, z0)

    def tr(l):
        return [1 - x for x in l]

    for op in ops:
        op._poly = Poly(op._poly.subs({z: 1/t, _x: -_x}), _x)
    slater2, cond2 = do_slater(tr(iq.bm), tr(iq.an), tr(iq.bq), tr(iq.ap),
                               t, 1/z0)

    slater1 = powdenest(slater1.subs(z, z0), polar=True)
    slater2 = powdenest(slater2.subs(t, 1/z0), polar=True)
    if not isinstance(cond2, bool):
        cond2 = cond2.subs(t, 1/z)

    m = meijerg(iq.an, iq.ap, iq.bm, iq.bq, z)
    if m.delta > 0 or \
       (m.delta == 0 and len(m.ap) == len(m.bq) and \
        (re(m.nu) < -1) is not False and polar_lift(z0) == polar_lift(1)):
        # The condition delta > 0 means that the convergence region is
        # connected. Any expression we find can be continued analytically
        # to the entire convergence region.
        # The conditions delta==0, p==q, re(nu) < -1 imply that G is continuous
        # on the positive reals, so the values at z=1 agree.
        if cond1 is not False:
            cond1 = True
        if cond2 is not False:
            cond2 = True

    if cond1 is True:
        slater1 = slater1.rewrite(rewrite or 'nonrep')
    else:
        slater1 = slater1.rewrite(rewrite or 'nonrepsmall')
    if cond2 is True:
        slater2 = slater2.rewrite(rewrite or 'nonrep')
    else:
        slater2 = slater2.rewrite(rewrite or 'nonrepsmall')

    if not isinstance(cond1, bool):
        cond1 = cond1.subs(z, z0)
    if not isinstance(cond2, bool):
        cond2 = cond2.subs(z, z0)

    def weight(expr, cond):
        from sympy import oo, zoo, nan
        if cond is True:
            c0 = 0
        elif cond is False:
            c0 = 1
        else:
            c0 = 2
        if expr.has(oo, zoo, -oo, nan):
            # XXX this actually should not happen, but consider
            # S('meijerg(((0, -1/2, 0, -1/2, 1/2), ()), ((0,),
            #   (-1/2, -1/2, -1/2, -1)), exp_polar(I*pi))/4')
            c0 = 3
        return (c0, expr.count(hyper), expr.count_ops())

    w1 = weight(slater1, cond1)
    w2 = weight(slater2, cond2)
    if min(w1, w2) <= (0, 1, oo):
        if w1 < w2:
            return slater1
        else:
            return slater2
    if max(w1[0], w2[0]) <= 1 and max(w1[1], w2[1]) <= 1:
        return Piecewise((slater1, cond1), (slater2, cond2),
                   (meijerg(iq_.an, iq_.ap, iq_.bm, iq_.bq, z0), True))

    # We couldn't find an expression without hypergeometric functions.
    # TODO it would be helpful to give conditions under which the integral
    #      is known to diverge.
    r = Piecewise((slater1, cond1), (slater2, cond2),
                   (meijerg(iq_.an, iq_.ap, iq_.bm, iq_.bq, z0), True))
    if r.has(hyper) and not allow_hyper:
        debug('  Could express using hypergeometric functions, '+
              'but not allowed.')
    if not r.has(hyper) or allow_hyper:
        return r

    return meijerg(iq_.an, iq_.ap, iq_.bm, iq_.bq, z0)


def hyperexpand(f, allow_hyper=False, rewrite='default'):
    """
    Expand hypergeometric functions. If allow_hyper is True, allow partial
    simplification (that is a result different from input,
    but still containing hypergeometric functions).

    Examples
    ========

    >>> from sympy.simplify.hyperexpand import hyperexpand
    >>> from sympy.functions import hyper
    >>> from sympy.abc import z
    >>> hyperexpand(hyper([], [], z))
    exp(z)

    Non-hyperegeometric parts of the expression and hypergeometric expressions
    that are not recognised are left unchanged:

    >>> hyperexpand(1 + hyper([1, 1, 1], [], z))
    hyper((1, 1, 1), (), z) + 1
    """
    from sympy.functions import hyper, meijerg
    from sympy import nan, zoo, oo
    f = sympify(f)

    def do_replace(ap, bq, z):
        r = _hyperexpand(IndexPair(ap, bq), z, rewrite=rewrite)
        if r is None:
            return hyper(ap, bq, z)
        else:
            return r

    def do_meijer(ap, bq, z):
        r = _meijergexpand(IndexQuadruple(ap[0], ap[1], bq[0], bq[1]), z,
                           allow_hyper, rewrite=rewrite)
        if not r.has(nan, zoo, oo, -oo):
            return r
    return f.replace(hyper, do_replace).replace(meijerg, do_meijer)

from sympy.polys.polytools import Poly
