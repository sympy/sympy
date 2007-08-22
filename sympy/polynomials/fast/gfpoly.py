"""Univariate polynomials with galois field coefficients."""

import random

import modint
import sparse_poly

def GFPolyFactory(p):
    """Create custom class for specific coefficient type."""
    coefficient_type = modint.ModularIntegerFactory(p)
    class newClass(sparse_poly.SparsePolynomial):
        coeff_type = coefficient_type
        zero = coeff_type(0)

        @staticmethod
        def from_int_dict(int_dict):
            """Alternative construction, through integers."""
            result_dict = {}
            for e, c in int_dict.iteritems():
                cc = coefficient_type(c)
                if cc:
                    result_dict[e] = cc
            return newClass(result_dict)

        def to_int_dict(self):
            """Returns the dictionaries of integer representators."""
            result_dict = {}
            for e, c in self.coeffs.iteritems():
                result_dict[e] = c.value
            return result_dict

        def to_sym_int_dict(self):
            """Returns the dictionaries of symmetric integer representators."""
            result_dict = {}
            for e, c in self.coeffs.iteritems():
                result_dict[e] = int(c)
            return result_dict

        @staticmethod
        def random(min_degree, max_degree, monic=True):
            """Generate random polynomial in given degree range."""
            degree = random.randrange(min_degree, max_degree + 1)
            p = coefficient_type.modulus
            result_dict = {}
            if monic:
                result_dict[degree] = coefficient_type(1)
                degree -= 1
            for e in xrange(0, degree + 1):
                c = coefficient_type(random.randrange(p))
                if c:
                    result_dict[e] = c
            return newClass(result_dict)
            
        def monic(self):
            if not self:
                return self.coeff_type(0), self
            leading_coeff = self[self.degree]
            return leading_coeff, self.scale(self.coeff_type(1)/leading_coeff)

    newClass.__name__ = "%sPoly" % coefficient_type.__name__
    return newClass


# Division algorithms:

def div(f, g):
    """Division with remainder."""
    q = f.__class__()
    r = f
    if not g:
        return q, r
    deg_diff = r.degree - g.degree
    while deg_diff >= 0:
        quot = f.__class__({deg_diff: r[r.degree]/g[g.degree]})
        q += quot
        r -= quot*g
        deg_diff = r.degree - g.degree
    return q, r

def gcd(f, g):
    """Euclidean algorithm."""
    while g:
        f, g = g, div(f,g)[1]
    return f.monic()[1]

def lcm(f, g):
    q, r = div(f*g, gcd(f,g))
    assert not r
    return q.monic()[1]

def xgcd(f, g):
    """Extended euclidean algorithm.

    Outputs the gcd, s and t, such that:
        h == s*f + t*g
    
    """
    one = f.coeff_type(1)
    p, q, r, s, t  = [], [], [], [], []
    pp, rr = f.monic()
    p.append(pp)
    r.append(rr)
    pp, rr = g.monic()
    p.append(pp)
    r.append(rr)
    s.append(f.__class__({0:(one/p[0])}))
    s.append(f.__class__())
    t.append(f.__class__())
    t.append(f.__class__({0:(one/p[1])}))
            
    while True:
        q.append(div(r[-2], r[-1])[0])
        pp, rr = (r[-2] - q[-1]*r[-1]).monic()
        if not rr:
            return r[-1], s[-1], t[-1]
        p.append(pp)
        r.append(rr)
        pp = one/pp
        s.append((s[-2] - q[-1]*s[-1]).scale(pp))
        t.append((t[-2] - q[-1]*t[-1]).scale(pp))

# Arithmetic modular a polynomial p:

def truncate(f, n):
    """The remainder from division by x**n."""
    result_dict = {}
    for e, c in f.coeffs.iteritems():
        if e < n:
            result_dict[e] = c
    return f.__class__(result_dict)


def pow_mod(f, n, p):
    """Repeated squaring."""
    assert isinstance(n, (int, long)) and n >= 0
    if n == 0:
        return f.__class__({0: f.__class__.coeff_type(1)})
    binary_n = []
    while n:
        if n % 2:
            binary_n.insert(0, 1)
            n = (n - 1) / 2
        else:
            binary_n.insert(0, 0)
            n /= 2
    result = div(f, p)[1]
    for k in binary_n[1:]:
        result *= result
        result = div(result, p)[1]
        if k:
            result *= f
            result = div(result, p)[1]
    return result


# Factorization:

def distinct_degree_factor(f):
    """Return a list of divisors.

    Each polynomial has only factors of a specific degree.
    """

    result = []
    coeff_type = f.__class__.coeff_type
    p = coeff_type.modulus
    x_poly = f.__class__({1: coeff_type(1)})
    one_poly = f.__class__({0: coeff_type(1)})
    h = x_poly
    while f != one_poly:
        h = pow_mod(h, p, f) # h <- h**p mod f
        g = gcd(h - x_poly, f)
        f, r = div(f, g)
        assert not r
        result.append(g)
        # Early abort:
        if f.degree < 2*(g.degree + 1):
            result.append(f)
            break
    return result

def equal_degree_split(f, degree):
    """Finds divisor of a result from distinct-degree factorization."""
    coeff_type = f.__class__.coeff_type
    one_poly = f.__class__({0: coeff_type(1)})
    a = f.random(1, f.degree - 1)
    g = gcd(f, a)
    if g != one_poly:
        return g
    b = pow_mod(a, (coeff_type.modulus**degree - 1)/2, f)
    g = gcd(b - one_poly, f)
    if g != one_poly and g != f:
        return g
    return None # Failure, try again with another random a.

def equal_degree_factor(f, degree):
    """Finds all divisors of a result from distinct-degree factorization."""
    if f.degree == degree:
        return [f]
    g = None
    while g is None:
        g = equal_degree_split(f, degree)
    q, r = div(f, g)
    assert not r
    return equal_degree_factor(g, degree) + equal_degree_factor(q, degree)

def factor(f):
    """Factorization of a univariate polynomial over a Galois field.

    Returns a list of the leading coefficient of f and the monic
    factors with their multiplicities.
    """
    p = f.__class__.coeff_type.modulus
    leading_coeff, f = f.monic()
    one_poly = f.__class__({0: f.__class__.coeff_type(1)})
    x_poly = f.__class__({1: f.__class__.coeff_type(1)})
    h = x_poly
    i = 0
    result = [leading_coeff]

    while f != one_poly:
        i += 1
        # One distinct-degree factorization step.
        h = pow_mod(h, p, f) # h <- h**p mod f
        g = gcd(h - x_poly, f)
        if g != one_poly:
            # Equal-degree factorization for degree i:
            g_factors = equal_degree_factor(g, i)
            # Now determine multiplicities of factors.
            for gg in g_factors:
                e = 0
                q, r = div(f, gg)
                while not r: # gg**e divides f
                    e += 1
                    f = q
                    q, r = div(f, gg)
                result.append((gg, e))
    return result

def factor_sqf(f):
    """Factorization of a univariate square-free polynomial over a Galois field.

    Returns a list of the leading coefficient and the monic factors of f.
    """

    one_poly = f.__class__({0: f.__class__.coeff_type(1)})
    leading_coeff, f = f.monic()
    result = [leading_coeff]
    for degree, divisor in enumerate(distinct_degree_factor(f)):
        if divisor == one_poly:
            continue
        result += equal_degree_factor(divisor, degree + 1)
    return result

