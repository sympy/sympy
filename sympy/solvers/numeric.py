# TODO: * calculate J numerically
#       * better exceptions for bad input
#       * solving overdetermined systems with Gauss-Newton algorithm

from __future__ import division

from sympy.matrices import Matrix
from sympy import Rational

# TODO: determine this from mpmath, currently assumes float precision
eps = 2**(-52) # = B/2*B**(-m)

def norm_p(x):
    """
    Return the p-norm of a vector.
    """
    return sum([abs(i)**p for i in x])**(1./p)

def maxnorm(x):
    """
    Return the oo-norm of a vector.
    """
    return max(map(abs, x))

def solve(A, b):
    """
    linear equation system Ax = b -> x
    """
    A = Matrix(A)
    b = Matrix(b)
    x = A.LUsolve(b)
    return x

def newton(f, x0, J, tol=eps, maxsteps=10, verbose=False, norm=maxnorm):
    """
    Find the root of a vector function numerically using Newton's method.

    f is a vector function representing a nonlinear equation system.
    x0 is the starting point close to the root.
    J is a function returning the jacobian matrix for a point.

    To be able to solve the system, it needs exactly as many variables as equations.

    Be careful with starting points. Using integers might lead to obscure
    fractions.
    """
    if isinstance(x0, (tuple, list)):
        x0 = Matrix(x0)
        assert x0.cols == 1, 'need a vector'
    fx = f(*x0)
    fxnorm = norm(fx)
    for k in xrange(maxsteps):
        if fxnorm <= tol:
            break
        # get direction of descent
        fxn = -fx
        Jx = J(*x0)
        s = solve(Jx, fxn)
        if verbose:
            print
            print 'step', str(k) + ':'
            print 'x0:', x0
            print 'fx:', fx
            print '||fx||:', fxnorm
            print 'Jx:'
            print Jx
            print 's:', s
        # damping step size TODO: better strategy (hard task)
        l = Rational(1)
        x1 = x0 + s
        while True:
            if x1 == x0:
                if verbose:
                    print "canceled, won't get more excact"
                return x0
            fx = f(*x1)
            newnorm = norm(fx)
            if verbose:
                print '||x1||:', newnorm
            if newnorm < fxnorm:
                # new x accepted
                fxnorm = newnorm
                x0 = x1
                break
            l /= 2
            x1 = x0 + l*s
            if verbose:
                print 'l:', l
                print 'x1:', x1
                print
    return x0
