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
from sympy import SYMPY_DEBUG

def add_formulae(formulae):
    """ Create our knowledge base.
        Leave this at the top for easy reference. """
    z = Dummy('z')
    a, b, c = symbols('a b c', cls=Dummy)
    def add(ap, bq, res):
        formulae.append(Formula(ap, bq, z, res, (a, b, c)))
    def addb(ap, bq, B, C, M):
        formulae.append(Formula(ap, bq, z, None, (a, b, c), B, C, M))

    from sympy.matrices import diag, Matrix

    # Luke, Y. L. (1969), The Special Functions and Their Approximations,
    # Volume 1, section 6.2

    from sympy import (exp, sqrt, cosh, log, asin, atan, I, lowergamma, cos,
                       atanh, besseli, gamma, erf, pi, sin, besselj)
    from sympy import polar_lift, exp_polar

    # 0F0
    add((), (), exp(z))

    # 1F0
    add((-a, ), (), (1+exp_polar(-I*pi)*z)**a)

    # 2F1
    addb((a, a - S.Half), (2*a,),
         Matrix([2**(2*a-1)*(1 + sqrt(1+exp_polar(-I*pi)*z))**(1-2*a),
                 2**(2*a-1)*(1 + sqrt(1+exp_polar(-I*pi)*z))**(-2*a)]),
         Matrix([[1, 0]]),
         Matrix([[(a-S.Half)*z/(1-z), (S.Half-a)*z/(1-z)],
                 [a/(1-z), a*(z-2)/(1-z)]]))
    addb((1, 1), (2,), # it is not particularly obvious, but this is the cts. branch
         Matrix([log(1 + exp_polar(-I*pi)*z), 1]), Matrix([[-1/z, 0]]),
         Matrix([[0, z/(z - 1)], [0, 0]]))
    # TODO branching
    addb((S.Half, 1), (S('3/2'),),
         Matrix([log((1 + sqrt(z))/(1 - sqrt(z)))/sqrt(z), 1]),
         Matrix([[S(1)/2, 0]]),
         Matrix([[-S(1)/2, 1/(1 - z)], [0, 0]]))
    addb((S.Half, S.Half), (S('3/2'),),
         Matrix([asin(sqrt(z))/sqrt(z), 1/sqrt(1 - z)]),
         Matrix([[1, 0]]),
         Matrix([[-S(1)/2, S(1)/2], [0, z/(1 - z)/2]]))
    # TODO branching
    addb((-a, S.Half - a), (S.Half,),
         Matrix([(1 + sqrt(z))**(2*a) + (1 - sqrt(z))**(2*a),
                 sqrt(z)*(1 + sqrt(z))**(2*a-1)
                 - sqrt(z)*(1 - sqrt(z))**(2*a-1)]),
         Matrix([[S.Half, 0]]),
         Matrix([[0, a], [z*(2*a-1)/2/(1-z), S.Half - z*(2*a-1)/(1-z)]]))

    # A. P. Prudnikov, Yu. A. Brychkov and O. I. Marichev (1990).
    # Integrals and Series: More Special Functions, Vol. 3,.
    # Gordon and Breach Science Publisher
    add([a, -a], [S.Half], cos(2*a*asin(sqrt(z))))
    # TODO branching
    addb([1, 1], [3*S.Half],
         Matrix([asin(sqrt(z))/sqrt(z*(1-z)), 1]), Matrix([[1, 0]]),
         Matrix([[(z - S.Half)/(1 - z), 1/(1 - z)/2], [0, 0]]))

    # 3F2
    # TODO branching
    addb([-S.Half, 1, 1], [S.Half, 2],
         Matrix([sqrt(z)*atanh(sqrt(z)), log(1 - z), 1]),
         Matrix([[-S(2)/3, -S(1)/(3*z), S(2)/3]]),
         Matrix([[S(1)/2, 0, z/(1 - z)/2],
                 [0, 0, z/(z - 1)],
                 [0, 0, 0]]))
    # actually the formula for 3/2 is much nicer ...
    # TODO branching
    addb([-S.Half, 1, 1], [2, 2],
         Matrix([sqrt(1 - z), log(sqrt(1 - z)/2 + S.Half), 1]),
         Matrix([[S(4)/9 - 16/(9*z), 4/(3*z), 16/(9*z)]]),
         Matrix([[z/2/(z - 1), 0, 0], [1/(2*(z - 1)), 0, S.Half], [0, 0, 0]]))

    # 1F1
    addb([1], [b], Matrix([z**(1 - b) * exp(z) * lowergamma(b - 1, z), 1]),
         Matrix([[b - 1, 0]]),Matrix([[1 - b + z, 1], [0, 0]]))
    addb([a], [2*a],
         Matrix([z**(S.Half - a)*exp(z/2)*besseli(a - S.Half, z/2)
                 * gamma(a + S.Half)/4**(S.Half - a),
                 z**(S.Half - a)*exp(z/2)*besseli(a + S.Half, z/2)
                 * gamma(a + S.Half)/4**(S.Half - a)]),
         Matrix([[1, 0]]),
         Matrix([[z/2, z/2], [z/2, (z/2 - 2*a)]]))
    add([-S.Half], [S.Half], exp(z) - sqrt(pi*z)*(-I)*erf(I*sqrt(z)))

    # 2F2
    addb([S.Half, a], [S(3)/2, a + 1],
         Matrix([a/(2*a - 1)*(-I)*sqrt(pi/z)*erf(I*sqrt(z)),
                 a/(2*a - 1)*(polar_lift(-1)*z)**(-a)*lowergamma(a, polar_lift(-1)*z), a/(2*a - 1)*exp(z)]),
         Matrix([[1, -1, 0]]),
         Matrix([[-S.Half, 0, 1], [0, -a, 1], [0, 0, z]]))

    # 0F1
    add((), (S.Half,), cosh(2*sqrt(z)))
    addb([], [b],
         Matrix([gamma(b)*z**((1-b)/2)*besseli(b-1, 2*sqrt(z)),
                 gamma(b)*z**(1 - b/2)*besseli(b  , 2*sqrt(z))]),
         Matrix([[1, 0]]), Matrix([[0, 1], [z, (1-b)]]))

    # 0F3
    x = 4*z**(S(1)/4)
    def fp(a,z): return besseli(a, x) + besselj(a, x)
    def fm(a,z): return besseli(a, x) - besselj(a, x)
    # TODO branching
    addb([], [S.Half, a, a+S.Half],
         Matrix([fp(2*a - 1, z), fm(2*a, z)*z**(S(1)/4),
                 fm(2*a - 1, z)*sqrt(z), fp(2*a, z)*z**(S(3)/4)])
           * 2**(-2*a)*gamma(2*a)*z**((1-2*a)/4),
         Matrix([[1, 0, 0, 0]]),
         Matrix([[0, 1, 0, 0],
                 [0, S(1)/2 - a, 1, 0],
                 [0, 0, S(1)/2, 1],
                 [z, 0, 0, 1 - a]]))
    x = 2*(-4*z)**(S(1)/4)
    addb([], [a, a + S.Half, 2*a],
         (2*sqrt(-z))**(1-2*a)*gamma(2*a)**2 *
         Matrix([besselj(2*a-1, x)*besseli(2*a-1, x),
                 x*(besseli(2*a, x)*besselj(2*a-1, x)
                    - besseli(2*a-1, x)*besselj(2*a, x)),
                 x**2*besseli(2*a, x)*besselj(2*a, x),
                 x**3*(besseli(2*a,x)*besselj(2*a-1,x)
                       + besseli(2*a-1, x)*besselj(2*a, x))]),
         Matrix([[1, 0, 0, 0]]),
         Matrix([[0, S(1)/4, 0, 0],
                 [0, (1-2*a)/2, -S(1)/2, 0],
                 [0, 0, 1-2*a, S(1)/4],
                 [-32*z, 0, 0, 1-a]]))

    # 1F2
    addb([a], [a - S.Half, 2*a],
         Matrix([z**(S.Half - a)*besseli(a-S.Half, sqrt(z))**2,
                 z**(1-a)*besseli(a-S.Half, sqrt(z))
                         *besseli(a-S(3)/2, sqrt(z)),
                 z**(S(3)/2-a)*besseli(a-S(3)/2, sqrt(z))**2]),
         Matrix([[-gamma(a + S.Half)**2/4**(S.Half - a),
                 2*gamma(a - S.Half)*gamma(a + S.Half)/4**(1 - a),
                 0]]),
         Matrix([[1 - 2*a, 1, 0], [z/2, S.Half - a, S.Half], [0, z, 0]]))
    addb([S.Half], [b, 2 - b],
         pi*(1-b)/sin(pi*b) *
         Matrix([besseli(1-b, sqrt(z))*besseli(b-1, sqrt(z)),
                 sqrt(z)*(besseli(-b, sqrt(z))*besseli(b-1, sqrt(z))
                          + besseli(1-b, sqrt(z))*besseli(b, sqrt(z))),
                 besseli(-b, sqrt(z))*besseli(b, sqrt(z))]),
         Matrix([[1, 0, 0]]),
         Matrix([[b-1, S(1)/2, 0],
                 [z, 0, z],
                 [0, S(1)/2, -b]]))

    # 2F3
    # XXX with this five-parameter formula is pretty slow with the current
    #     Formula.find_instantiations (creates 2!*3!*3**(2+3) ~ 3000
    #     instantiations ... But it's not too bad.
    addb([a, a + S.Half], [2*a, b, 2*a - b + 1],
         gamma(b)*gamma(2*a - b + 1) * (sqrt(z)/2)**(1-2*a) *
         Matrix([besseli(b-1, sqrt(z))*besseli(2*a-b, sqrt(z)),
                 sqrt(z)*besseli(b, sqrt(z))*besseli(2*a-b, sqrt(z)),
                 sqrt(z)*besseli(b-1, sqrt(z))*besseli(2*a-b+1, sqrt(z)),
                 besseli(b, sqrt(z))*besseli(2*a-b+1, sqrt(z))]),
         Matrix([[1, 0, 0, 0]]),
         Matrix([[0, S(1)/2, S(1)/2, 0],
                 [z/2, 1-b, 0, z/2],
                 [z/2, 0, b-2*a, z/2],
                 [0, S(1)/2, S(1)/2, -2*a]]))

def add_meijerg_formulae(formulae):
    z = Dummy('z')
    a, b, c = map(Dummy, 'abc')
    rho = Dummy('rho')
    def add(an, ap, bm, bq, B, C, M, matcher):
        formulae.append(MeijerFormula(an, ap, bm, bq, z, [a, b, c, rho],
                                      B, C, M, matcher))

    from sympy import Matrix
    from sympy import gamma, uppergamma, exp

    def detect_uppergamma(iq):
        x = iq.an[0]
        y, z = iq.bm
        swapped = False
        if Mod1(x) == Mod1(y):
            swapped = True
            (y, z) = (z, y)
        if Mod1(x) != Mod1(z) or x > z:
            return None
        l = [y, x]
        if swapped:
            l = [x, y]
        return {rho: y, a: x-y}, IndexQuadruple([x], [], l, [])
    add([a+rho], [], [rho, a+rho], [],
        Matrix([gamma(1 - a)*z**rho*exp(z)*uppergamma(a, z),
                gamma(1-a)*z**(a+rho)]),
        Matrix([[1, 0]]),
        Matrix([[rho+z, -1], [0, a+rho]]),
        detect_uppergamma)


def make_simp(z):
    """ Create a function that simplifies rational functions in `z`. """
    def simp(expr):
        """ Efficiently simplify the rational function `expr`. """
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


class Mod1(object):
    """
    Represent an expression 'mod 1'.

    Beware: __eq__ and the hash are NOT compatible. (by design)
    This means that m1 == m2 does not imply hash(m1) == hash(m2).
    Code that creates Mod1 objects (like compute_buckets below) should be
    careful only to produce one instance of Mod1 for each class.
    """
    # TODO this should be backported to any implementation of a Mod object
    #      (c/f issue 2490)

    def __new__(cls, r):
        if r.is_Rational and not r.free_symbols:
            return r - r.p//r.q
        res = object.__new__(cls)
        res.expr = r
        return res

    def __repr__(self):
        return str(self.expr) + ' % 1'

    #Needed to allow adding Mod1 objects to a dict in Python 3
    def __hash__(self):
        return super(Mod1, self).__hash__()

    def __eq__(self, other):
        from sympy import simplify
        if not isinstance(other, Mod1):
            return False
        if simplify(self.expr - other.expr).is_integer is True:
            return True
        return False

    def __ne__(self, other):
        return not self == other

class IndexPair(object):
    """ Holds a pair of indices, and methods to compute their invariants. """

    def __init__(self, ap, bq):
        from sympy import expand, Tuple
        self.ap = Tuple(*[expand(x) for x in sympify(ap)])
        self.bq = Tuple(*[expand(x) for x in sympify(bq)])

    @property
    def sizes(self):
        return (len(self.ap), len(self.bq))

    def __str__(self):
        return 'IndexPair(%s, %s)' % (self.ap, self.bq)

    def compute_buckets(self, oabuckets=None, obbuckets=None):
        """
        Partition parameters `ap`, `bq` into buckets, that is return two dicts
        abuckets, bbuckets such that every key in [ab]buckets is a rational in
        range [0, 1) and the corresponding items are items of ap/bq congruent to
        the key mod 1.

        If oabuckets, obbuckets is specified, try to use the same Mod1 objects
        for parameters where possible.

        >>> from sympy.simplify.hyperexpand import IndexPair
        >>> from sympy import S
        >>> ap = (S(1)/2, S(1)/3, S(-1)/2, -2)
        >>> bq = (1, 2)
        >>> IndexPair(ap, bq).compute_buckets()
        ({0: (-2,), 1/3: (1/3,), 1/2: (1/2, -1/2)}, {0: (1, 2)})
        """
        # TODO this should probably be cached somewhere
        abuckets = {}
        bbuckets = {}

        oaparametric = []
        obparametric = []
        if oabuckets is not None:
            for parametric, buckets in [(oaparametric, oabuckets),
                                        (obparametric, obbuckets)]:
                parametric += filter(lambda x: isinstance(x, Mod1),
                                     buckets.keys())

        for params, bucket, oparametric in [(self.ap, abuckets, oaparametric),
                                            (self.bq, bbuckets, obparametric)]:
            parametric = []
            for p in params:
                res = Mod1(p)
                if isinstance(res, Mod1):
                    parametric.append(p)
                    continue
                if res in bucket:
                    bucket[res] += (p,)
                else:
                    bucket[res] = (p,)
            while parametric:
                p0 = parametric[0]
                p0mod1 = Mod1(p0)
                if oparametric.count(p0mod1):
                    i = oparametric.index(p0mod1)
                    p0mod1 = oparametric.pop(i)
                bucket[p0mod1] = (p0,)
                pos = []
                for po in parametric[1:]:
                    if Mod1(po) == p0mod1:
                        bucket[p0mod1] += (po,)
                    else:
                        pos.append(po)
                parametric = pos

        return abuckets, bbuckets

    def build_invariants(self):
        """
        Compute the invariant vector of (`ap`, `bq`), that is:
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
            if not any(isinstance(x[0], Mod1) for x in bucket):
                bucket.sort(key=lambda x: x[0])
            bucket = tuple(map(lambda x: (x[0], len(x[1])), bucket))
            return bucket

        return (gamma, tr(abuckets), tr(bbuckets))

    def difficulty(self, ip):
        """ Estimate how many steps it takes to reach `ip` from self.
            Return -1 if impossible. """
        oabuckets, obbuckets = self.compute_buckets()
        abuckets, bbuckets = ip.compute_buckets(oabuckets, obbuckets)

        gt0 = lambda x: (x > 0) is True
        if S(0) in abuckets and (not S(0) in oabuckets or
             len(filter(gt0, abuckets[S(0)])) != len(filter(gt0, oabuckets[S(0)]))):
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
        def tr(l): return Tuple(*[expand(x) for x in sympify(l)])
        self.an = tr(an)
        self.ap = tr(ap)
        self.bm = tr(bm)
        self.bq = tr(bq)

    def compute_buckets(self):
        """
        Compute buckets for the fours sets of parameters.
        We guarantee that any two equal Mod1 objects returned are actually the
        same, and that the buckets are sorted by real part (an and bq
        descendending, bm and ap ascending).

        >>> from sympy.simplify.hyperexpand import IndexQuadruple
        >>> from sympy.abc import y
        >>> from sympy import S
        >>> IndexQuadruple([1, 3, 2, S(3)/2], [1 + y, y, 2, y + 3], [2], [y]).compute_buckets()
        ({0: [3, 2, 1], 1/2: [3/2]}, {0: [2], y + 1 % 1: [y, y + 1, y + 3]}, {0: [2]}, {y + 1 % 1: [y]})

        """
        mod1s = []
        pan, pap, pbm, pbq = {}, {}, {}, {}
        for dic, lis in [(pan, self.an), (pap, self.ap), (pbm, self.bm),
                           (pbq, self.bq)]:
            for x in lis:
                m = Mod1(x)
                if mod1s.count(m):
                    i = mod1s.index(m)
                    m = mod1s[i]
                else:
                    mod1s.append(m)
                dic.setdefault(m, []).append(x)

        for dic, flip in [(pan, True), (pap, False), (pbm, False), (pbq, True)]:
            l = dic.items()
            dic.clear()
            for m, items in l:
                x0 = items[0]
                items.sort(key=lambda x: x-x0)
                if flip:
                    items.reverse()
                dic[m] = items

        return pan, pap, pbm, pbq

    @property
    def signature(self):
        return (len(self.an), len(self.ap), len(self.bm), len(self.bq))

    def __str__(self):
        return 'IndexQuadruple(%s, %s, %s, %s)' % (self.an, self.ap,
                                                   self.bm, self.bq)

# Dummy generator
x = Dummy('x')

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

        afactors = map(lambda a: x + a, self.indices.ap)
        bfactors = map(lambda b: x + b - 1, self.indices.bq)
        expr = x*Mul(*bfactors) - self.z*Mul(*afactors)
        poly = Poly(expr, x)

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
        ap = Tuple(*map(expand, sympify(ap)))
        bq = Tuple(*map(expand, sympify(bq)))
        z  = sympify(z)
        res = sympify(res)
        symbols = filter(lambda x: ap.has(x) or bq.has(x), sympify(symbols))

        self.z  = z
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
                    if m != (a,) or not c.is_Rational:
                        raise NotImplementedError('?')
                    l = ilcm(l, c.q)

                    if not p.has(*others):
                        isolating.append((i, c.q, c.p))
                lcms[a] = l
                i += 1
            if len(isolating) == 0:
                raise NotImplementedError('parameter is not isolated')
            isolating.sort(key=lambda x:x[1])
            isolating.sort(key=lambda x:-x[2])
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
        `ip.ap`, `ip.bq`. Return the instantiated formulae as a list.
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
        Decide if `self` is a suitable origin.

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
        >>> Formula((S(1)/2, 1), (2, -S(2)/3, S(3)/2), None, None, []).is_suitable()
        True
        """
        from sympy import oo, zoo
        if len(self.symbols) > 0:
            return None
        for a in self.indices.ap:
            for b in self.indices.bq:
                if (a-b).is_integer and not a < b:
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
        Given the suitable parameters `ip.ap`, `ip.bq`, try to find an origin
        in our knowledge base.

        >>> from sympy.simplify.hyperexpand import FormulaCollection, IndexPair
        >>> f = FormulaCollection()
        >>> f.lookup_origin(IndexPair((), ())).closed_form
        exp(_z)
        >>> f.lookup_origin(IndexPair([1], ())).closed_form
        1/(_z*exp_polar(-I*pi) + 1)

        >>> from sympy import S
        >>> f.lookup_origin(IndexPair([S('1/4'), S('3/4 + 4')], [S.Half])).closed_form
        1/(2*(sqrt(_z) + 1)**(17/2)) + 1/(2*(-sqrt(_z) + 1)**(17/2))
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
        possible.sort(key=lambda x:x[0])
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
        an = Tuple(*map(expand, sympify(an)))
        ap = Tuple(*map(expand, sympify(ap)))
        bm = Tuple(*map(expand, sympify(bm)))
        bq = Tuple(*map(expand, sympify(bq)))
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
        formulae = []
        add_meijerg_formulae(formulae)
        self.formulae = {}
        for formula in formulae:
            self.formulae.setdefault(formula.indices.signature, []).append(formula)

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
        Apply `self` to the object `obj`, where the generator is given by `op`.

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
        diffs  = [obj]
        for c in coeffs[1:]:
            diffs.append(op(diffs[-1]))
        r = coeffs[0]*diffs[0]
        for c, d in zip(coeffs[1:], diffs[1:]):
            r += c*d
        return r

class MultOperator(Operator):
    """ Simply multiply by a "constant" """

    def __init__(self, p):
        self._poly = Poly(p, x)

class ShiftA(Operator):
    """ Increment an upper index. """

    def __init__(self, ai):
        ai = sympify(ai)
        if ai == 0:
            raise ValueError('Cannot increment zero upper index.')
        self._poly = Poly(x/ai + 1, x)

    def __str__(self):
        return '<Increment upper %s.>' % (1/self._poly.all_coeffs()[0])

class ShiftB(Operator):
    """ Decrement a lower index. """

    def __init__(self, bi):
        bi = sympify(bi)
        if bi == 1:
            raise ValueError('Cannot decrement unit lower index.')
        self._poly = Poly(x/(bi - 1) + 1, x)

    def __str__(self):
        return '<Decrement lower %s.>' % (1/self._poly.all_coeffs()[0] + 1)

class UnShiftA(Operator):
    """ Decrement an upper index. """

    def __init__(self, ap, bq, i, z):
        """ Note: i counts from zero! """
        ap, bq, i = map(sympify, [ap, bq, i])

        self._ap = ap
        self._bq = bq
        self._i  = i

        ap = list(ap)
        bq = list(bq)
        ai = ap.pop(i) - 1

        if ai == 0:
            raise ValueError('Cannot decrement unit upper index.')

        m = Poly(z*ai, x)
        for a in ap:
            m *= Poly(x + a, x)
        #print m

        A = Dummy('A')
        D = Poly(ai*A - ai, A)
        n = 1*D
        for b in bq:
            n *= (D + b - 1)
        #print n

        b0 = -n.nth(0)
        if b0 == 0:
            raise ValueError('Cannot decrement upper index: ' \
                               'cancels with lower')
        #print b0

        n = Poly(Poly(n.all_coeffs()[:-1], A).as_expr().subs(A, x/ai + 1), x)

        self._poly = Poly((n-m)/b0, x)

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
        self._i  = i

        ap = list(ap)
        bq = list(bq)
        bi = bq.pop(i) + 1

        if bi == 0:
            raise ValueError('Cannot increment -1 lower index.')

        m = Poly(x*(bi-1), x)
        for b in bq:
            m *= Poly(x + b - 1, x)
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

        n = Poly(Poly(n.all_coeffs()[:-1], B).as_expr().subs(B, x/(bi-1) + 1), x)
        #print n

        self._poly = Poly((m-n)/b0, x)

    def __str__(self):
        return '<Increment lower index #%s of %s, %s.>' % (self._i,
                                                        self._ap, self._bq)

class MeijerShiftA(Operator):
    """ Increment an upper b index. """

    def __init__(self, bi):
        bi = sympify(bi)
        self._poly = Poly(bi - x, x)

    def __str__(self):
        return '<Increment upper b=%s.>' % (self._poly.all_coeffs()[1])

class MeijerShiftB(Operator):
    """ Decrement an upper a index. """

    def __init__(self, bi):
        bi = sympify(bi)
        self._poly = Poly(1 - bi + x, x)

    def __str__(self):
        return '<Decrement upper a=%s.>' % (1 - self._poly.all_coeffs()[1])

class MeijerShiftC(Operator):
    """ Increment a lower b index. """

    def __init__(self, bi):
        bi = sympify(bi)
        self._poly = Poly(-bi + x, x)

    def __str__(self):
        return '<Increment lower b=%s.>' % (-self._poly.all_coeffs()[1])

class MeijerShiftD(Operator):
    """ Decrement a lower a index. """

    def __init__(self, bi):
        bi = sympify(bi)
        self._poly = Poly(bi - 1 - x, x)

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
        self._i  = i

        an = list(an)
        ap = list(ap)
        bm = list(bm)
        bq = list(bq)
        bi = bm.pop(i) - 1

        m = Poly(1, x)
        for b in bm:
            m *= Poly(b - x, x)
        for b in bq:
            m *= Poly(x - b, x)
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

        n = Poly(Poly(n.all_coeffs()[:-1], A).as_expr().subs(A, bi - x), x)
        #print n

        self._poly = Poly((m-n)/b0, x)

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
        self._i  = i

        an = list(an)
        ap = list(ap)
        bm = list(bm)
        bq = list(bq)
        ai = an.pop(i) + 1

        m = Poly(z, x)
        for a in an:
            m *= Poly(1 - a + x, x)
        for a in ap:
            m *= Poly(a - 1 - x, x)
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

        n = Poly(Poly(n.all_coeffs()[:-1], B).as_expr().subs(B, 1 - ai + x), x)
        #print n

        self._poly = Poly((m-n)/b0, x)

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
        self._i  = i

        an = list(an)
        ap = list(ap)
        bm = list(bm)
        bq = list(bq)
        bi = bq.pop(i) - 1

        m = Poly(1, x)
        for b in bm:
            m *= Poly(b - x, x)
        for b in bq:
            m *= Poly(x - b, x)
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

        n = Poly(Poly(n.all_coeffs()[:-1], C).as_expr().subs(C, x - bi), x)
        #print n

        self._poly = Poly((m-n)/b0, x)

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
        self._i  = i

        an = list(an)
        ap = list(ap)
        bm = list(bm)
        bq = list(bq)
        ai = ap.pop(i) + 1

        m = Poly(z, x)
        for a in an:
            m *= Poly(1 - a + x, x)
        for a in ap:
            m *= Poly(a - 1 - x, x)
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

        n = Poly(Poly(n.all_coeffs()[:-1], B).as_expr().subs(B, ai - 1 - x), x)
        #print n

        self._poly = Poly((m-n)/b0, x)

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
            p *= (x + bj + k)/(bj + k)

        self._poly = Poly(p, x)
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
            p *= (sign*x + a + k)

        self._poly = Poly(p, x)
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
    """ Order reduction algorithm common to both Hypergeometric and Meijer G """
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
    Given the hypergeometric parameters `ip.ap`, `ip.bq`, find a sequence of operators
    to reduces order as much as possible.

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
    (IndexPair((2,), (3,)), [<Reduce order by cancelling upper 4 with lower 3.>])
    """
    nap, nbq, operators = _reduce_order(ip.ap, ip.bq, ReduceOrder, lambda x: x)

    return IndexPair(Tuple(*nap), Tuple(*nbq)), operators

def reduce_order_meijer(iq):
    """
    Given the Meijer G function parameters, `iq.am`, `iq.ap`, `iq.bm`,
    `iq.bq`, find a sequence of operators that reduces order as much as possible.

    Return niq, [operators].

    Examples
    ========

    >>> from sympy.simplify.hyperexpand import reduce_order_meijer, IndexQuadruple
    >>> reduce_order_meijer(IndexQuadruple([3, 4], [5, 6], [3, 4], [1, 2]))[0]
    IndexQuadruple((4, 3), (5, 6), (3, 4), (2, 1))
    >>> reduce_order_meijer(IndexQuadruple([3, 4], [5, 6], [3, 4], [1, 8]))[0]
    IndexQuadruple((3,), (5, 6), (3, 4), (1,))
    >>> reduce_order_meijer(IndexQuadruple([3, 4], [5, 6], [7, 5], [1, 5]))[0]
    IndexQuadruple((3,), (), (), (1,))
    >>> reduce_order_meijer(IndexQuadruple([3, 4], [5, 6], [7, 5], [5, 3]))[0]
    IndexQuadruple((), (), (), ())
    """

    nan, nbq, ops1 = _reduce_order(iq.an, iq.bq, ReduceOrder.meijer_plus, lambda x: -x)
    nbm, nap, ops2 = _reduce_order(iq.bm, iq.ap, ReduceOrder.meijer_minus, lambda x: x)

    return IndexQuadruple(Tuple(*nan), Tuple(*nap), Tuple(*nbm), Tuple(*nbq)), \
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
    Apply the list of operators `ops` to object `obj`, substituting `op` for the
    generator.
    """
    res = obj
    for o in reversed(ops):
        res = o.apply(res, op)
    return res

def devise_plan(ip, nip, z):
    """
    Devise a plan (consisting of shift and un-shift operators) to be applied
    to the hypergeometric function (`nip.ap`, `nip.bq`) to yield
    (`ip.ap`, `ip.bq`).
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
    >>> devise_plan(IndexPair((1, S.Half), ()), IndexPair((2, S('3/2')), ()), z)
    [<Decrement upper index #0 of [2, 1/2], [].>, <Decrement upper index #0 of [3/2, 2], [].>]

    A slightly more complicated plan:

    >>> devise_plan(IndexPair((1, 3), ()), IndexPair((2, 2), ()), z)
    [<Increment upper 2.>, <Decrement upper index #0 of [2, 2], [].>]

    Another more complicated plan: (note that the ap have to be shifted first!)

    >>> devise_plan(IndexPair((1, -1), (2,)), IndexPair((3, -2), (4,)), z)
    [<Decrement lower 3.>, <Decrement lower 4.>, <Decrement upper index #1 of [-1, 2], [4].>, <Decrement upper index #1 of [-1, 3], [4].>, <Increment upper -2.>]
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

    for r in set(abuckets.keys() + bbuckets.keys()):
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
            raise ValueError('%s not reachable from %s' % ((ap, bq), (nap, nbq)))

        al = sorted(list(al))
        nal = sorted(list(nal))
        bk = sorted(list(bk))
        nbk = sorted(list(nbk))

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
            amax  = al[-1]

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
    from sympy import oo, factorial, rf
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
    for n in xrange(-a):
        fac *= z
        fac /= n + 1
        for a in ip.ap: fac *= a + n
        for b in ip.bq: fac /= b + n
        res += fac
    return res

collection = None
def _hyperexpand(ip, z, ops0=[], z0=Dummy('z0'), premult=1, prem=0):
    """
    Try to find an expression for the hypergeometric function
    `ip.ap`, `ip.bq`.

    The result is expressed in terms of a dummy variable z0. Then it
    is multiplied by premult. Then ops0 is applied.
    premult must be a*z**prem for some a independent of z.
    """
    from sympy.simplify import powdenest, simplify, polarify
    z = polarify(z, subs=False)

    # TODO
    # The following would be possible:
    # 1) Partial simplification (i.e. return a simpler hypergeometric function,
    #    even if we cannot express it in terms of named special functions).
    # 2) PFD Duplication (see Kelly Roach's paper)
    # 3) If the coefficients are a rational function of n (numerator parameters
    #    k, a1, ..., an, denominator parameters a1+k1, a2+k2, ..., an+kn, where
    #    k, k1, ..., kn are integers) then result can be expressed using Lerch
    #    transcendent. Under certain conditions, this simplifies to polylogs
    #    or even zeta functions. C/f Kelly Roach's paper.

    global collection
    if collection is None:
        collection = FormulaCollection()

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
        return simplify(p).subs(z0, z)

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

    # Now try to find a formula
    f = collection.lookup_origin(nip)

    if f is None:
        debug('  Could not find an origin.')
        # There is nothing we can do.
        return None

    # We need to find the operators that convert f into (nap, nbq).
    ops += devise_plan(nip, f.indices, z0)

    # Now carry out the plan.
    C = apply_operators(f.C.subs(f.z, z0), ops,
                        make_derivative_operator(f.M.subs(f.z, z0), z0))
    from sympy import eye
    C = apply_operators(C, ops0,
                        make_derivative_operator(f.M.subs(f.z, z0)
                                                 + prem*eye(f.M.shape[0]), z0))

    if premult == 1:
        C = C.applyfunc(make_simp(z0))
    r = C*f.B.subs(f.z, z0)*premult
    r = r[0].subs(z0, z) + p

    # This will simpliy things like sqrt(-z**2) to i*z.
    # It would be wrong under certain choices of branch, but all results we
    # return are under an "implicit suitable choice of branch" anyway.
    return powdenest(r, polar=True) # This unpolarifies!

def devise_plan_meijer(fro, to, z):
    """
    Find a sequence of operators to convert index quadruple `fro` into
    index quadruple `to`. It is assumed that fro and to have the same
    signatures, and that in fact any corresponding pair of parameters differs
    by integers, and a direct path is possible. I.e. if there are parameters
       a1 b1 c1  and a2 b2 c2
    it is assumed that a1 can be shifted to a2, etc.
    The only thing this routine determines is the order of shifts to apply,
    nothing clever will be tried.
    It is also assumed that fro is suitable.

    >>> from sympy.simplify.hyperexpand import devise_plan_meijer, IndexQuadruple
    >>> from sympy.abc import z

    Empty plan:

    >>> devise_plan_meijer(IndexQuadruple([1], [2], [3], [4]), IndexQuadruple([1], [2], [3], [4]), z)
    []

    Very simple plans:

    >>> devise_plan_meijer(IndexQuadruple([0], [], [], []), IndexQuadruple([1], [], [], []), z)
    [<Increment upper a index #0 of [0], [], [], [].>]
    >>> devise_plan_meijer(IndexQuadruple([0], [], [], []), IndexQuadruple([-1], [], [], []), z)
    [<Decrement upper a=0.>]
    >>> devise_plan_meijer(IndexQuadruple([], [1], [], []), IndexQuadruple([], [2], [], []), z)
    [<Increment lower a index #0 of [], [1], [], [].>]

    Slightly more complicated plans:

    >>> devise_plan_meijer(IndexQuadruple([0], [], [], []), IndexQuadruple([2], [], [], []), z)
    [<Increment upper a index #0 of [1], [], [], [].>, <Increment upper a index #0 of [0], [], [], [].>]
    >>> devise_plan_meijer(IndexQuadruple([0], [], [0], []), IndexQuadruple([-1], [], [1], []), z)
    [<Increment upper b=0.>, <Decrement upper a=0.>]

    Order matters:

    >>> devise_plan_meijer(IndexQuadruple([0], [], [0], []), IndexQuadruple([1], [], [1], []), z)
    [<Increment upper a index #0 of [0], [], [1], [].>, <Increment upper b=0.>]
    """
    # TODO for now, we use the following simple heuristic: inverse-shift
    #      when possible, shift otherwise. Give up if we cannot make progress.
    def try_shift(f, t, shifter, diff, counter):
        """ Try to apply `shifter` in order to bring some element in `f` nearer
            to its counterpart in `to`. `diff` is +/- 1 and determines the
            effect of `shifter`. Counter is a list of elements blocking the
            shift.
            Return an operator if change was possible, else None.
        """
        for idx, (a, b) in enumerate(zip(f, t)):
            if (
                (a-b).is_integer and (b-a)/diff > 0 and
                all(a != x for x in counter)
            ):
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

meijercollection = None
def _meijergexpand(iq, z0, allow_hyper=False):
    """
    Try to find an expression for the Meijer G function specified
    by the IndexQuadruple `iq`. If `allow_hyper` is True, then returning
    an expression in terms of hypergeometric functions is allowed.

    Currently this just does slater's theorem.
    """
    from sympy import hyper, Piecewise, meijerg, powdenest
    global meijercollection
    if meijercollection is None:
        meijercollection = MeijerFormulaCollection()

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
    f = meijercollection.lookup_origin(iq)
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
    # 2) Paired Index Theorems
    # 3) PFD Duplication
    #    (See Kelly Roach's paper for (2) and (3).)
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

    def do_slater(an, bm, ap, bq, z):
        from sympy import gamma, residue, factorial, rf, expand_func, polar_lift
        iq = IndexQuadruple(an, bm, ap, bq)
        _, pbm, pap, _ = iq.compute_buckets()
        if not can_do(pbm, pap):
            return S(0), False

        res = S(0)
        for m in pbm:
            if len(pbm[m]) == 1:
                bh = pbm[m][0]
                fac = 1
                bo = list(bm)
                bo.remove(bh)
                for bj in bo: fac *= gamma(bj - bh)
                for aj in an: fac *= gamma(1 + bh - aj)
                for bj in bq: fac /= gamma(1 + bh - bj)
                for aj in ap: fac /= gamma(aj - bh)
                nap = [1 + bh - a for a in list(an) + list(ap)]
                nbq = [1 + bh - b for b in list(bo) + list(bq)]

                k = polar_lift(S(-1)**(len(ap) - len(bm)))
                harg = k*z
                # NOTE even though k "is" +-1, this has to be t/k instead of
                #      t*k ... we are using polar numbers for consistency!
                premult = (t/k)**bh
                hyp = _hyperexpand(IndexPair(nap, nbq), harg, ops,
                                   t, premult, bh)
                if hyp is None:
                    hyp = apply_operators(premult*hyper(nap, nbq, t), ops,
                                          lambda f: t*f.diff(t)).subs(t, harg)
                res += fac * hyp
            else:
                b_ = pbm[m][0]
                ki = [bi - b_ for bi in pbm[m][1:]]
                u = len(ki)
                li = [ai - b_ for ai in pap[m][0:u+1]]
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
                harg = k*z
                premult = (t/k)**au
                nap = [1 + au - a for a in list(an) + list(ap)] + [1]
                nbq = [1 + au - b for b in list(bm) + list(bq)]

                hyp = _hyperexpand(IndexPair(nap, nbq), harg, ops,
                                   t, premult, au)
                if hyp is None:
                    hyp = apply_operators(premult*hyper(nap, nbq, t), ops,
                                          lambda f: t*f.diff(t)).subs(t, harg)

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

        cond = len(an) + len(ap) < len(bm) + len(bq)
        if len(an) + len(ap) == len(bm) + len(bq):
            cond = abs(z) < 1
        return res, cond

    t = Dummy('t')
    slater1, cond1 = do_slater(iq.an, iq.bm, iq.ap, iq.bq, z)

    def tr(l): return [1 - x for x in l]
    for op in ops:
        op._poly = Poly(op._poly.subs({z: 1/t, x: -x}), x)
    slater2, cond2 = do_slater(tr(iq.bm), tr(iq.an), tr(iq.bq), tr(iq.ap),
                               t)

    slater1 = powdenest(slater1.subs(z, z0), polar=True)
    slater2 = powdenest(slater2.subs(t, 1/z0), polar=True)
    if not isinstance(cond2, bool):
        cond2 = cond2.subs(t, 1/z0)

    if meijerg(iq.an, iq.ap, iq.bm, iq.bq, z).delta > 0:
        # The above condition means that the convergence region is connected.
        # Any expression we find can be continued analytically to the entire
        # convergence region.
        if cond1 is not False:
            cond1 = True
        if cond2 is not False:
            cond2 = True

    if not isinstance(cond1, bool): cond1 = cond1.subs(z, z0)
    if not isinstance(cond2, bool): cond2 = cond2.subs(z, z0)

    if cond1 is True and cond2 is True and not slater1.has(hyper) and \
       not slater2.has(hyper):
        # If both are possible and free of unevaluated terms, return the
        # simpler one.
        from sympy import count_ops
        if count_ops(slater1) > count_ops(slater2):
            return slater2
        else:
            return slater1

    if cond1 is True and not slater1.has(hyper):
        return slater1
    if cond2 is True and not slater2.has(hyper):
        return slater2

    # We couldn't find an expression without hypergeometric functions.
    # TODO it would be helpful to give conditions under which the integral
    #      is known to diverge.
    r =  Piecewise((slater1, cond1), (slater2, cond2),
                   (meijerg(iq_.an, iq_.ap, iq_.bm, iq_.bq, z0), True))
    if r.has(hyper) and not allow_hyper:
        debug('  Could express using hypergeometric functions, but not allowed.')
    if not r.has(hyper) or allow_hyper:
        return r

    return meijerg(iq_.an, iq_.ap, iq_.bm, iq_.bq, z0)

def hyperexpand(f, allow_hyper=False):
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
        r = _hyperexpand(IndexPair(ap, bq), z)
        if r is None:
            return hyper(ap, bq, z)
        else:
            return r
    def do_meijer(ap, bq, z):
        r = _meijergexpand(IndexQuadruple(ap[0], ap[1], bq[0], bq[1]), z,
                           allow_hyper)
        if not r.has(nan, zoo, oo, -oo):
            return r
    return f.replace(hyper, do_replace).replace(meijerg, do_meijer)

from sympy.polys.polytools import Poly
