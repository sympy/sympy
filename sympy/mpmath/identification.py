"""
Implements the PSLQ algorithm for integer relation detection,
and derivative algorithms for constant recognition.
"""

from mptypes import (mpf, eps, mp)
from functions import (log, exp, sqrt)

from libmpf import to_fixed, from_man_exp, MODE
from libelefun import mpf_sqrt
from libelefun import sqrt_fixed as _sqrt_fixed

# round to nearest integer (can be done more elegantly...)
def round_fixed(x, prec):
    return ((x + (1<<(prec-1))) >> prec) << prec

# XXX. this is needed because the Python sqrt_fixed currently chokes on
# x not ~= 1
if MODE == 'python':
    def sqrt_fixed(x, prec):
        return to_fixed(mpf_sqrt(from_man_exp(x, -prec, prec), prec), prec)
else:
    sqrt_fixed = _sqrt_fixed

def pslq(x, eps=None):
    """
    Given a vector of real numbers x = [x1, x2, ..., xn], pslq(x) uses the
    PSLQ algorithm to find a list of integers [c1, c2, ..., cn] such that
    c1*x1 + c2*x2 + ... + cn*xn = 0 approximately.


    This is a fairly direct translation to Python of the pseudocode given by
    David Bailey, "The PSLQ Integer Relation Algorithm":
    http://www.cecm.sfu.ca/organics/papers/bailey/paper/html/node3.html

    The stopping criteria are NOT yet properly implemented.

    Note: now using fixed-point arithmetic for a ~7x speedup compared
    to the original, pure-mpf version.
    """
    n = len(x)
    assert n >= 1
    prec = mp.prec
    assert prec >= 53
    target = prec // max(2,n)
    if target < 30:
        if target < 5:
            print "Warning: precision for PSLQ may be too low"
        target = int(prec * 0.75)
    if eps is None:
        eps = mpf(2)**(-target)
    extra = 60
    prec += extra
    eps = to_fixed(eps._mpf_, prec)
    x = [None] + [to_fixed(mpf(xk)._mpf_, prec) for xk in x]
    g = sqrt_fixed((4<<prec)//3, prec)
    A = {}
    B = {}
    H = {}
    # Initialization
    # step 1
    for i in range(1, n+1):
        for j in range(1, n+1):
            A[i,j] = B[i,j] = (i==j) << prec
            H[i,j] = 0
    # step 2
    s = [None] + [0] * n
    for k in range(1, n+1):
        t = 0
        for j in range(k, n+1):
            t += (x[j]**2 >> prec)
        s[k] = sqrt_fixed(t, prec)
    t = s[1]
    y = x[:]
    for k in range(1, n+1):
        y[k] = (x[k] << prec) // t
        s[k] = (s[k] << prec) // t
    # step 3
    for i in range(1, n+1):
        for j in range(i+1, n):
            H[i,j] = 0
        if i <= n-1:
            H[i,i] = (s[i+1] << prec) // s[i]
        for j in range(1, i):
            H[i,j] = ((-y[i]*y[j])<<prec)//(s[j]*s[j+1])
    # step 4
    for i in range(2, n+1):
        for j in range(i-1, 0, -1):
            #t = floor(H[i,j]/H[j,j] + 0.5)
            t = round_fixed((H[i,j] << prec)//H[j,j], prec)
            y[j] = y[j] + (t*y[i] >> prec)
            for k in range(1, j+1):
                H[i,k] = H[i,k] - (t*H[j,k] >> prec)
            for k in range(1, n+1):
                A[i,k] = A[i,k] - (t*A[j,k] >> prec)
                B[k,j] = B[k,j] + (t*B[k,i] >> prec)
    # Main algorithm
    for REP in range(100):
        # step 1
        m = -1
        szmax = -1
        for i in range(1, n):
            h = H[i,i]
            sz = (sqrt_fixed((4<<prec)//3, prec)**i * abs(h)) >> (prec*(i-1))
            if sz > szmax:
                m = i
                szmax = sz
        # step 2
        y[m], y[m+1] = y[m+1], y[m]
        tmp = {}
        for i in range(1,n+1): H[m,i], H[m+1,i] = H[m+1,i], H[m,i]
        for i in range(1,n+1): A[m,i], A[m+1,i] = A[m+1,i], A[m,i]
        for i in range(1,n+1): B[i,m], B[i,m+1] = B[i,m+1], B[i,m]
        # step 3
        if m <= n - 2:
            t0 = sqrt_fixed((H[m,m]**2 + H[m,m+1]**2)>>prec, prec)
            # XXX: this could be spurious, due to fixed-point arithmetic
            if not t0:
                break
            t1 = (H[m,m] << prec) // t0
            t2 = (H[m,m+1] << prec) // t0
            for i in range(m, n+1):
                t3 = H[i,m]
                t4 = H[i,m+1]
                H[i,m] = (t1*t3+t2*t4) >> prec
                H[i,m+1] = (-t2*t3+t1*t4) >> prec
        # step 4
        for i in range(m+1, n+1):
            for j in range(min(i-1, m+1), 0, -1):
                try:
                    t = round_fixed((H[i,j] << prec)//H[j,j], prec)
                # XXX
                except ZeroDivisionError:
                    break
                y[j] = y[j] + ((t*y[i]) >> prec)
                for k in range(1, j+1):
                    H[i,k] = H[i,k] - (t*H[j,k] >> prec)
                for k in range(1, n+1):
                    A[i,k] = A[i,k] - (t*A[j,k] >> prec)
                    B[k,j] = B[k,j] + (t*B[k,i] >> prec)
        for i in range(1, n+1):
            if abs(y[i]) < eps:
                vec = [int(round_fixed(B[j,i], prec) >> prec) for j in range(1,n+1)]
                if max(abs(v) for v in vec) < 10**6:
                    return vec
    return None


def findpoly(x, n=1):
    """Find an integer polynomial P of degree at most n such that
    x is an approximate root of P."""
    if x == 0:
        return [1, 0]
    xs = [mpf(1)]
    for i in range(1,n+1):
        xs.append(x**i)
        a = pslq(xs)
        if a is not None:
            return a[::-1]

def fracgcd(p, q):
    x, y = p, q
    while y:
        x, y = y, x % y
    if x != 1:
        p //= x
        q //= x
    if q == 1:
        return p
    return p, q

def pslqstring(r, constants):
    q = r[0]
    r = r[1:]
    s = []
    for i in range(len(r)):
        p = r[i]
        if p:
            z = fracgcd(-p,q)
            cs = constants[i][1]
            if cs == '1':
                cs = ''
            else:
                cs = '*' + cs
            if isinstance(z, (int, long)):
                if z > 0: term = str(z) + cs
                else:     term = ("(%s)" % z) + cs
            else:
                term = ("(%s/%s)" % z) + cs
            s.append(term)
    s = ' + '.join(s)
    if '+' in s or '*' in s:
        s = '(' + s + ')'
    return s or '0'

def prodstring(r, constants):
    q = r[0]
    r = r[1:]
    num = []
    den = []
    for i in range(len(r)):
        p = r[i]
        if p:
            z = fracgcd(-p,q)
            cs = constants[i][1]
            if isinstance(z, (int, long)):
                if abs(z) == 1: t = cs
                else:           t = '%s**%s' % (cs, abs(z))
                ([num,den][z<0]).append(t)
            else:
                t = '%s**(%s/%s)' % (cs, abs(z[0]), z[1])
                ([num,den][z[0]<0]).append(t)
    num = '*'.join(num)
    den = '*'.join(den)
    if num and den: return "(%s)/(%s)" % (num, den)
    if num: return num
    if den: return "1/(%s)" % den

def quadraticstring(t,a,b,c):
    if c < 0:
        a,b,c = -a,-b,-c
    u1 = (-b+sqrt(b**2-4*a*c))/(2*c)
    u2 = (-b-sqrt(b**2-4*a*c))/(2*c)
    if abs(u1-t) < abs(u2-t):
        if b:  s = '((%s+sqrt(%s))/%s)' % (-b,b**2-4*a*c,2*c)
        else:  s = '(sqrt(%s)/%s)' % (-4*a*c,2*c)
    else:
        if b:  s = '((%s-sqrt(%s))/%s)' % (-b,b**2-4*a*c,2*c)
        else:  s = '(-sqrt(%s)/%s)' % (-4*a*c,2*c)
    return s

# Transformation y = f(x,c), with inverse function x = f(y,c)
# The third entry indicates whether the transformation is
# redundant when c = 1
transforms = [
  (lambda x,c: x*c, '$y/$c', 0),
  (lambda x,c: x/c, '$c*$y', 1),
  (lambda x,c: c/x, '$c/$y', 0),
  (lambda x,c: (x*c)**2, 'sqrt($y)/$c', 0),
  (lambda x,c: (x/c)**2, '$c*sqrt($y)', 1),
  (lambda x,c: (c/x)**2, '$c/sqrt($y)', 0),
  (lambda x,c: c*x**2, 'sqrt($y)/sqrt($c)', 1),
  (lambda x,c: x**2/c, 'sqrt($c)*sqrt($y)', 1),
  (lambda x,c: c/x**2, 'sqrt($c)/sqrt($y)', 1),
  (lambda x,c: sqrt(x*c), '$y**2/$c', 0),
  (lambda x,c: sqrt(x/c), '$c*$y**2', 1),
  (lambda x,c: sqrt(c/x), '$c/$y**2', 0),
  (lambda x,c: c*sqrt(x), '$y**2/$c**2', 1),
  (lambda x,c: sqrt(x)/c, '$c**2*$y**2', 1),
  (lambda x,c: c/sqrt(x), '$c**2/$y**2', 1),
  (lambda x,c: exp(x*c), 'log($y)/$c', 0),
  (lambda x,c: exp(x/c), '$c*log($y)', 1),
  (lambda x,c: exp(c/x), '$c/log($y)', 0),
  (lambda x,c: c*exp(x), 'log($y/$c)', 1),
  (lambda x,c: exp(x)/c, 'log($c*$y)', 1),
  (lambda x,c: c/exp(x), 'log($c/$y)', 0),
  (lambda x,c: log(x*c), 'exp($y)/$c', 0),
  (lambda x,c: log(x/c), '$c*exp($y)', 1),
  (lambda x,c: log(c/x), '$c/exp($y)', 0),
  (lambda x,c: c*log(x), 'exp($y/$c)', 1),
  (lambda x,c: log(x)/c, 'exp($c*$y)', 1),
  (lambda x,c: c/log(x), 'exp($c/$y)', 0),
]

def identify(x, constants=[], full=False, maxcoeff=1000, tolerance=None,
    quadratics=True, verbose=False):
    """
    This function attempts to find a symbolic expression for the given
    quantity x. It can identify simple algebraic numbers, as well as
    simple combinations of the given list of base constants, and
    exponentials or logarithms thereof.

    The base constants should be given as strings representing mpmath
    expressions.

    If a match is found, a mathematical formula is returned as a string.
    With full=True, a list of matching formulas is returned.

    In order not to produce spurious approximations, high precision
    should be used; preferrably 50 digits or more.

    Examples:

        >>> mp.dps = 15
        >>> identify(0.22222222222222222)
        '(2/9)'

        >>> mp.dps = 50
        >>> identify(3*pi + 4*sqrt(2), ['pi','sqrt(2)'])
        ...
        (3*pi + 4*sqrt(2))

    Further example identifications that should work (many redundant
    results may be found if run with full=True):

        mp.dps = 50
        base = ['sqrt(2)','pi','log(2)']
        identify(0.25, base)
        identify(3*pi + 2*sqrt(2) + 5*log(2)/7, base)
        identify(exp(pi+2), base)
        identify(1/(3+sqrt(2)), base)
        identify(sqrt(2)/(3*pi+4), base)
        identify(5**(mpf(1)/3)*pi*log(2)**2, base)
    """

    solutions = []

    def addsolution(s):
        if verbose: print "Found: ", s
        solutions.append(s)

    x = mpf(x)

    # Further along, x will be assumed positive
    if x == 0:
        if full: return ['0']
        else:    return '0'

    if x < 0:
        sol = identify(-x, constants, full, maxcoeff, tolerance, quadratics, verbose)
        if sol is None:
            return sol
        if full:
            return ["-(%s)"%s for s in sol]
        else:
            return "-(%s)" % sol

    sols = []
    if tolerance:
        weps = mpf(tolerance)
    else:
        weps = eps**0.7

    if isinstance(constants, dict):
        constants = [(mpf(v), name) for (name, v) in constants.items()]
    else:
        import sympy.mpmath
        constants = [(eval(p, sympy.mpmath.__dict__), p) for p in constants]

    # We always want to find at least rational terms
    if 1 not in [value for (name, value) in constants]:
        constants = [(mpf(1), '1')] + constants

    # PSLQ with simple algebraic and functional transformations
    for ft, ftn, red in transforms:
        for c, cn in constants:
            if red and cn == '1':
                continue
            t = ft(x,c)
            # Prevent exponential transforms from wreaking havoc
            if abs(t) > maxcoeff**2 or abs(t) < weps:
                continue
            # Linear combination of base constants
            r = pslq([t] + [a[0] for a in constants], weps)
            s = None
            if r is not None and max(abs(uw) for uw in r) <= maxcoeff and r[0]:
                s = pslqstring(r, constants)
            # Quadratic algebraic numbers
            elif quadratics:
                q = pslq([mpf(1), t, t**2], weps)
                if q is not None and len(q) == 3 and q[2]:
                    aa, bb, cc = q
                    if max(abs(aa),abs(bb),abs(cc)) <= maxcoeff:
                        s = quadraticstring(t,aa,bb,cc)
            if s:
                if cn == '1' and ('/$c' in ftn):
                    s = ftn.replace('$y', s).replace('/$c', '')
                else:
                    s = ftn.replace('$y', s).replace('$c', cn)
                addsolution(s)
                if not full: return solutions[0]

            if verbose:
                print "."

    # Check for a direct multiplicative formula
    if x != 1:
        # Allow fractional powers of fractions
        ilogs = [2,3,5,7]
        # Watch out for existing fractional powers of fractions
        logs = []
        for a, s in constants:
            if not sum(bool(findpoly(log(a)/log(i),1)) for i in ilogs):
                logs.append((log(a), s))
        logs = [(log(i),str(i)) for i in ilogs] + logs
        r = pslq([log(x)] + [a[0] for a in logs], weps)
        if r is not None and max(abs(uw) for uw in r) <= maxcoeff and r[0]:
            addsolution(prodstring(r, logs))
            if not full: return solutions[0]

    if full:
        return sorted(solutions, key=len)
    else:
        return None
