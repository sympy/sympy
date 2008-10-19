# -*- encoding: utf-8 -*-

from mptypes import convert_lossless, extraprec, eps
from calculus import diff, diffc
from functions import sqrt, sign

##############
# 1D-SOLVERS #
##############

class Secant:
    """
    1d-solver generating pairs of approximative root and error.

    Needs starting points x0 and x1 close to the root.
    x1 defaults to x0 + 0.25.

    Pro:
    * converges fast
    Contra:
    * converges slowly for multiple roots
    """
    maxsteps = 30

    def __init__(self, f, x0, **kwargs):
        if len(x0) == 1:
            self.x0 = x0[0]
            self.x1 = self.x0 + 0.25
        elif len(x0) == 2:
            self.x0 = x0[0]
            self.x1 = x0[1]
        else:
            raise ValueError('expected 1 or 2 starting points, got %i' * len(x0))
        self.f = f

    def __iter__(self):
        f = self.f
        x0 = self.x0
        x1 = self.x1
        f0 = f(x0)
        while True:
            f1 = f(x1)
            l = x1 - x0
            if not l:
                break
            s = (f1 - f0) / l
            if not s:
                break
            x0, x1 = x1, x1 - f1/s
            f0 = f1
            yield x1, abs(l)

class MNewton:
    """
    1d-solver generating pairs of approximative root and error.

    Needs starting point x0 close to the root.
    Uses modified Newton's method that converges fast regardless of the
    multiplicity of the root.

    Pro:
    * converges fast for multiple roots
    Contra:
    * needs first and second derivative of f
    * 3 function evaluations per iteration
    """
    maxsteps = 20

    def __init__(self, f, x0, **kwargs):
        if not len(x0) == 1:
            raise ValueError('expected 1 starting point, got %i' * len(x0))
        self.x0 = x0[0]
        self.f = f
        if not 'df' in kwargs:
            def df(x):
                return diff(f, x)
        else:
            df = kwargs['df']
        self.df = df
        if not 'd2f' in kwargs:
            def d2f(x):
                return diff(df, x)
        else:
            d2f = kwargs['df']
        self.d2f = d2f

    def __iter__(self):
        x = self.x0
        f = self.f
        df = self.df
        d2f = self.d2f
        while True:
            prevx = x
            fx = f(x)
            if fx == 0:
                break
            dfx = df(x)
            d2fx = d2f(x)
            # x = x - F(x)/F'(x) with F(x) = f(x)/f'(x)
            x -= fx / (dfx - fx * d2fx / dfx)
            error = abs(x - prevx)
            yield x, error

class Halley:
    """
    1d-solver generating pairs of approximative root and error.

    Needs a starting point x0 close to the root.
    Uses Halley's method with cubic convergance rate.

    Pro:
    * converges even faster the Newton's method
    * useful when computing with *many* digits
    Contra:
    * needs first and second derivative of f
    * 3 function evaluations per iteration
    * converges slowly for multiple roots
    """

    maxsteps = 20

    def __init__(self, f, x0, **kwargs):
        if not len(x0) == 1:
            raise ValueError('expected 1 starting point, got %i' * len(x0))
        self.x0 = x0[0]
        self.f = f
        if not 'df' in kwargs:
            def df(x):
                return diff(f, x)
        else:
            df = kwargs['df']
        self.df = df
        if not 'd2f' in kwargs:
            def d2f(x):
                return diff(df, x)
        else:
            d2f = kwargs['df']
        self.d2f = d2f

    def __iter__(self):
        x = self.x0
        f = self.f
        df = self.df
        d2f = self.d2f
        while True:
            prevx = x
            fx = f(x)
            dfx = df(x)
            d2fx = d2f(x)
            x -=  2*fx*dfx / (2*dfx**2 - fx*d2fx)
            error = abs(x - prevx)
            yield x, error

class Muller:
    """
    1d-solver generating pairs of approximative root and error.

    Needs starting points x0, x1 and x2 close to the root.
    x1 defaults to x0 + 0.25; x2 to x1 + 0.25.
    Uses Muller's method that converges towards complex roots.

    Pro:
    * converges fast (somewhat faster than secant)
    * can find complex roots
    Contra:
    * converges slowly for multiple roots
    * may have complex values for real starting points and real roots

    http://en.wikipedia.org/wiki/M%C3%BCller%27s_method
    """
    maxsteps = 30

    def __init__(self, f, x0, **kwargs):
        if len(x0) == 1:
            self.x0 = x0[0]
            self.x1 = self.x0 + 0.25
            self.x2 = self.x1 + 0.25
        elif len(x0) == 2:
            self.x0 = x0[0]
            self.x1 = x0[1]
            self.x2 = self.x1 + 0.25
        elif len(x0) == 3:
            self.x0 = x0[0]
            self.x1 = x0[1]
            self.x2 = x0[2]
        else:
            raise ValueError('expected 1, 2 or 3 starting points, got %i'
                             % len(x0))
        self.f = f
        self.verbose = kwargs['verbose']

    def __iter__(self):
        f = self.f
        x0 = self.x0
        x1 = self.x1
        x2 = self.x2
        fx0 = f(x0)
        fx1 = f(x1)
        fx2 = f(x2)
        while True:
            # TODO: maybe refactoring with function for divided differences
            # calculate divided diffferences
            fx2x1 = (fx1 - fx2) / (x1 - x2)
            fx2x0 = (fx0 - fx2) / (x0 - x2)
            fx1x0 = (fx0 - fx1) / (x0 - x1)
            w = fx2x1 + fx2x0 - fx1x0
            fx2x1x0 = (fx1x0 - fx2x1) / (x0 - x2)
            if w == 0 and fx2x1x0 == 0:
                if self.verbose:
                    print 'canceled with'
                    print 'x0 =', x0, ', x1 =', x1, 'and x2 =', x2
                break
            x0 = x1
            fx0 = fx1
            x1 = x2
            fx1 = fx2
            # denominator should be as large as possible => choose sign
            r = sqrt(w**2 - 4*fx2*fx2x1x0)
            if abs(w - r) > abs(w + r):
                r = -r
            x2 -= 2*fx2 / (w + r)
            fx2 = f(x2)
            error = abs(x2 - x1)
            yield x2, error

# TODO: consider raising a ValueError when there's no sign change in a and b
class Bisection:
    """
    1d-solver generating pairs of approximative root and error.

    Uses bisection method to find a root of f in [a, b].
    Might fail for multiple roots (needs sign change).

    Pro:
    * robust and reliable
    Contra:
    * converges slowly
    * needs sign change
    """
    maxsteps = 100

    def __init__(self, f, x0, **kwargs):
        if len(x0) != 2:
            raise ValueError('expected interval of 2 points, got %i' * len(x0))
        self.f = f
        self.a = x0[0]
        self.b = x0[1]

    def __iter__(self):
        f = self.f
        a = self.a
        b = self.b
        l = b - a
        while True:
            m = 0.5 * (a + b) # TODO: ldexp?
            if f(m) * f(b) < 0:
                a = m
            else:
                b = m
            l /= 2
            yield (a + b)/2, abs(l)

def _getm(method):
    """
    Return a function to calculate m for Illinois-like methods.
    """
    if method == 'illinois':
        def getm(fz, fb):
            return 0.5
    elif method == 'pegasus':
        def getm(fz, fb):
            return fb/(fb + fz)
    elif method == 'anderson':
        def getm(fz, fb):
            m = 1 - fz/fb
            if m > 0:
                return m
            else:
                return 0.5
    else:
        raise ValueError, "method '%s' not recognized" % method
    return getm

class Illinois:
    """
    1d-solver generating pairs of approximative root and error.

    Uses Illinois method or similair to find a root of f in [a, b].
    Might fail for multiple roots (needs sign change).
    Combines bisect with secant (improved regula falsi).

    The only difference between the methods is the scaling factor m, which is
    used to ensure convergence (you can choose one using the 'method' keyword):
    Illinois method ('illinois'):       m = 0.5
    Pegasus method ('pegasus'):         m = fb/(fb + fz)
    Anderson-Björk method ('anderson'): m = 1 - fz/fb if positive else 0.5

    Pro:
    * converges very fast
    Contra:
    * has problems with multiple roots
    * needs sign change
    """
    maxsteps = 30

    def __init__(self, f, x0, **kwargs):
        if len(x0) != 2:
            raise ValueError('expected interval of 2 points, got %i' * len(x0))
        self.a = x0[0]
        self.b = x0[1]
        self.f = f
        self.tol = kwargs['tol']
        self.verbose = kwargs['verbose']
        self.method = kwargs.get('method', 'illinois')
        self.getm = _getm(self.method)
        if self.verbose:
            print 'using %s method' % self.method

    def __iter__(self):
        method = self.method
        f = self.f
        a = self.a
        b = self.b
        fa = f(a)
        fb = f(b)
        m = None
        while True:
            l = b - a
            s = (fb - fa) / l
            z = a - fa/s
            fz = f(z)
            if abs(fz) < self.tol:
                # TODO: better condition (when f is very flat)
                if self.verbose:
                    print 'canceled with z =', z
                yield z, l
                break
            if fz * fb < 0: # root in [z, b]
                a = b
                fa = fb
                b = z
                fb = fz
            else: # root in [a, z]
                m = getm(fz, fb)
                b = z
                fb = fz
                fa = m*fa # scale down to ensure convergence
            if self.verbose and m and not method == 'illinois':
                print 'm:', m
            yield (a + b)/2, abs(l)

def Pegasus(*args, **kwargs):
    """
    1d-solver generating pairs of approximative root and error.

    Uses Pegasus method to find a root of f in [a, b].
    Wrapper for illinois to use method='pegasus'.
    """
    kwargs['method'] = 'pegasus'
    return Illinois(*args, **kwargs)

def Anderson(*args, **kwargs):
    u"""
    1d-solver generating pairs of approximative root and error.

    Uses Anderson-Björk method to find a root of f in [a, b].
    Wrapper for illinois to use method='pegasus'.
    """
    kwargs['method'] = 'anderson'
    return Illinois(*args, **kwargs)

# TODO: check whether it's possible to combine it with Illinois stuff
class Ridder:
    """
    1d-solver generating pairs of approximative root and error.

    Ridders' method to find a root of f in [a, b].
    Is told to perform as well as Brent's method while being simpler.

    Pro:
    * very fast
    * simpler than Brent's method
    Contra:
    * two function evaluations per step
    * has problems with multiple roots
    * needs sign change

    http://en.wikipedia.org/wiki/Ridders%27_method
    """
    maxsteps=30

    def __init__(self, f, x0, **kwargs):
        self.f = f
        if len(x0) != 2:
            raise ValueError('expected interval of 2 points, got %i' * len(x0))
        self.x1 = x0[0]
        self.x2 = x0[1]
        self.verbose = kwargs['verbose']
        self.tol = kwargs['tol']

    def __iter__(self):
        f = self.f
        x1 = self.x1
        fx1 = f(x1)
        x2 = self.x2
        fx2 = f(x2)
        while True:
            x3 = 0.5*(x1 + x2)
            fx3 = f(x3)
            x4 = x3 + (x3 - x1) * sign(fx1 - fx2) * fx3 / sqrt(fx3**2 - fx1*fx2)
            fx4 = f(x4)
            if abs(fx4) < self.tol:
                # TODO: better condition (when f is very flat)
                if self.verbose:
                    print 'canceled with f(x4) =', fx4
                yield x4, abs(x1 - x2)
                break
            if fx4 * fx2 < 0: # root in [x4, x2]
                x1 = x4
                fx1 = fx4
            else: # root in [x1, x4]
                x2 = x4
                fx2 = fx4
            error = abs(x1 - x2)
            yield (x1 + x2)/2, error

class ANewton:
    """
    EXPERIMENTAL 1d-solver generating pairs of approximative root and error.

    Uses Newton's method modified to use Steffensens method when convergence is
    slow. (I. e. for multiple roots.)
    """
    maxsteps = 20

    def __init__(self, f, x0, **kwargs):
        if not len(x0) == 1:
            raise ValueError('expected 1 starting point, got %i' * len(x0))
        self.x0 = x0[0]
        self.f = f
        if not 'df' in kwargs:
            def df(x):
                return diff(f, x)
        else:
            df = kwargs['df']
        self.df = df
        def phi(x):
            return x - f(x) / df(x)
        self.phi = phi
        self.verbose = kwargs['verbose']

    def __iter__(self):
        x0 = self.x0
        f = self.f
        df = self.df
        phi = self.phi
        error = 0
        counter = 0
        while True:
            prevx = x0
            try:
                x0 = phi(x0)
            except ZeroDivisionError:
                if self.verbose:
                    'ZeroDivisionError: canceled with x =', x0
                break
            preverror = error
            error = abs(prevx - x0)
            # TODO: decide not to use convergence acceleration
            if error and abs(error - preverror) / error < 1:
                if self.verbose:
                    print 'converging slowly'
                counter += 1
            if counter >= 3:
                # accelerate convergence
                phi = steffensen(phi)
                counter = 0
                if self.verbose:
                    print 'accelerating convergence'
            yield x0, error

# TODO: add Brent

#############
# UTILITIES #
#############

str2solver = {'secant':Secant,'mnewton':MNewton, 'halley':Halley,
              'muller':Muller, 'bisect':Bisection, 'illinois':Illinois,
              'pegasus':Pegasus, 'anderson':Anderson, 'ridder':Ridder,
              'anewton':ANewton}

@extraprec(20)
def findroot(f, x0, solver=Secant, tol=None, verbose=False, verify=True,
             force_type=convert_lossless, **kwargs):
    """
    Find a root of f using x0 as starting point or interval.

    If not abs(f(root)) < tol an exception is raised.

    Arguments:
    f : one dimensional function
    x0 : starting point, several starting points or interval (depends on solver)
    tol : the returned solution has an error smaller than this
    verbose : print additional information for each iteration if true
    verify : verify the solution and raise a ValueError if abs(f(x)) > tol
    force_type : use specified type constructor on starting points
    solver : a generator for f and x0 returning approximative solution and error
    maxsteps : after how many steps the solver will cancel
    df : first derivative of f (used by some solvers)
    d2f : second derivative of f (used by some solvers)

    solver has to be callable with (f, x0, **kwargs) and return an generator
    yielding pairs of approximative solution and estimated error.
    You can use the following string aliases:
    'secant', 'mnewton', 'halley', 'muller', 'illinois', 'pegasus', 'anderson',
    'ridder', 'anewton', 'bisect'
    See mpmath.optimization for their documentation.
    """
    # initialize arguments
    if not force_type:
        force_type = lambda x: x
    elif not tol and (force_type == float or force_type == complex):
        tol = 2**(-52) # TODO: consider a less strict value
    kwargs['verbose'] = verbose
    if 'd1f' in kwargs:
        kwargs['df'] = kwargs['d1f']
    if tol is None:
        tol = eps * 2**10
    kwargs['tol'] = tol
    if isinstance(x0, (list, tuple)):
        x0 = [force_type(x) for x in x0]
    else:
        x0 = [force_type(x0)]
    if isinstance(solver, str):
        try:
            solver = str2solver[solver]
        except KeyError:
            raise ValueError('could not recognize solver')
    # use solver
    iterations = solver(f, x0, **kwargs)
    if 'maxsteps' in kwargs:
        maxsteps = kwargs['maxsteps']
    else:
        maxsteps = iterations.maxsteps
    i = 0
    for x, error in iterations:
        if verbose:
            print 'x:', x
            print 'error:', error
        i += 1
        if error < tol or i >= maxsteps:
            break
    if verify and abs(f(x))**2 > tol: # TODO: better condition?
        raise ValueError('Could not find root within given tolerance. '
                         '(%g > %g)\n'
                         'Try another starting point or tweak arguments.'
                         % (abs(f(x)), tol))
    return x

def multiplicity(f, root, tol=eps, maxsteps=10, **kwargs):
    """
    Return the multiplicity of a given root of f.

    Internally, numerical derivatives are used. This is very inefficient for
    higher order derviatives. You can be specify the n-th derivative using the
    dnf keyword.
    """
    kwargs['d0f'] = f
    for i in xrange(maxsteps):
        dfstr = 'd' + str(i) + 'f'
        if dfstr in kwargs:
            df = kwargs[dfstr]
        else:
            df = lambda x: diffc(f, x, i)
        if not abs(df(root)) < tol:
            break
    return i

def steffensen(f):
    """
    linear convergent function -> quadratic convergent function

    Steffensen's method for quadratic convergence of a linear converging
    sequence.
    Don not use it for higher rates of convergence.
    It may even work for divergent sequences.

    Definition:
    F(x) = (x*f(f(x)) - f(x)**2) / (f(f(x)) - 2*f(x) + x)
    """
    def F(x):
        fx = f(x)
        ffx = f(fx)
        return (x*ffx - fx**2) / (ffx - 2*fx + x)
    return F
