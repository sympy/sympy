""" Hyperbolastic functions """
from sympy.functions.elementary.hyperbolic import asinh, tanh
from sympy.functions.elementary.exponential import exp
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.core.sympify import sympify
from sympy.core.function import Function, ArgumentIndexError

class h1(Function):
    r""" 
    The hyperbolastic rate equation of type I, denoted H1, is given by:

        dP(x)/dx = P(x)/M*(M-P(x))(delta + theta/(1+x**2)**(1/2))

    where x is any real number and P(x) is the population size at x.
    The parameter M represents carrying capacity, and parameters delta
    and theta  jointly represent growth rate. The parameter theta gives the
    distance from a symmetric sigmoidal curve.

    Examples
    ========

    >>> from sympy import diff, symbols, h1
    >>> a ,b, c, x, y, z = symbols('M delta theta x0 P0 x')

    >>> h1(a, b, c, x, y, z)
    M*P0/(P0 + (M - P0)*exp(-delta*x + delta*x0 - theta*asinh(x) - theta*asinh(x0)))

    >>> h1(a, b, 0, 1, 8, z)
    M*exp(delta*x)/((M - 8)*exp(delta)/8 + exp(delta*x))

    >>> diff(h1(a, b, 0, x, y, z), z).simplify()
    M*P0*delta*(M - P0)*exp(delta*(x + x0))/(P0*exp(delta*x) + (M - P0)*exp(delta*x0))**2

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Hyperbolastic_functions
    .. [2] https://tbiomed.biomedcentral.com/articles/10.1186/1742-4682-2-14

    """

    def fdiff(self, argindex=1):
        if argindex == 6:
            M, delta, theta, x0, P0 , x= self.args
            p_x = h1(M, delta, theta, x0, P0 , x)
            print(p_x)
            return (p_x/M*(M - p_x)*(delta + theta/(sqrt(1 + x**2)))).simplify()
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls ,M, delta, theta, x0, P0 , x):
        alpha = (M - P0)/P0*exp(delta*x0 - theta*asinh(x0))
        return M/(1 + alpha*exp(-delta*x - theta*asinh(x))).simplify()

class h2(Function):
    r""" 
    The hyperbolastic rate equation of type II, denoted H2, is given by:

    dP(x)/dx = alpha*delta*gamma*P(x)**2*x**(gamma - 1)/M*tanh(M - P(x)/alpha*p(x))

    with initial condition P(t0) = P0 and γ > 0, where tanh stands for hyperbolic
    tangent function, M is the carrying capacity, and beta and gamma are parameters.
    As in the H1 model, parameter beta has to be positive for increasing growth
    curves with an asymptote at M and is negative only for decay profiles.
    We refer to the growth rate model above as the hyperbolastic differential
    equation of type II.


    Examples
    ========

    >>> from sympy import diff, h2, symbols
    >>> a ,b, c, x, y, z = symbols('M delta gamma x0 P0 x')

    >>> h2(a, b, c, x, y, z)
    M*P0/(P0 + (M - P0)*asinh(exp(-delta*gamma*x))*asinh(exp(-delta*gamma*x0)))

    >>> h2(a, 0, c, x, y, z)
    M*P0/(P0 + (M - P0)*log(1 + sqrt(2))**2)

    >>> diff(h2(a, a, a, a, y, z), z)
    M**3*P0*(M - P0)*exp(-M**2*x)*asinh(exp(-M**3))/(sqrt(1 + exp(-2*M**2*x))*(P0 + (M - P0)*asinh(exp(-M**3))*asinh(exp(-M**2*x)))**2)

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Hyperbolastic_functions
    .. [2] https://tbiomed.biomedcentral.com/articles/10.1186/1742-4682-2-14

    """

    def fdiff(self, argindex=1):
        if argindex == 6:
            M, delta, gama, x0, P0 , x= self.args
            p_x = h2(M, delta, gama, x0, P0 , x)
            alpha = (M - P0)/P0*asinh(exp(-delta*x0*gama))
            return (p_x**2*x**(gama-1)*gama*delta*alpha)/M*tanh((M - p_x) \
                                                    /alpha*p_x).simplify()
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls ,M, delta, gama, x0, P0 , x):
            if sympify(gama).is_real and gama <= 0 :
                return ValueError("gamma must be greater than 0")
            else:
                alpha = (M - P0)/P0*asinh(exp(-delta*x0*gama))
                return M/(1 + alpha*asinh(exp(-delta*x*gama))).simplify()


class h3(Function):
    r"""
    The hyperbolastic rate equation of type III, denoted H3, is given by:

        dP(t)/dt = (M - P(t))*(delta*gamma*t**(gamma - 1) + theta/(1 + (theta*t)**2)**(1/2))

    with initial condition P(t0) = P0 where M is the carrying capacity t is the
    time and delta, gamma and theta are parameters. We refer to above model as
    the hyperbolastic ordinary differential equation of type III.

    Examples
    ========

    >>> from sympy import diff, h3, symbols
    >>> a ,b, c, d, x, y, z = symbols('M delta gamma theta t0 P0 t')

    >>> h3(a, b, c, d, x, y, z)
    M - (M - P0)*exp(-delta*t**gamma - asinh(t*theta))*exp(delta*t0**gamma + asinh(t0*theta))

    >>> h3(a, 0, c, d, x, y, z)
    M - (M - P0)*exp(-asinh(t*theta))*exp(asinh(t0*theta))

    >>> diff(h3(a, b, c, 0, x, y, z), z).simplify()
    delta*gamma*t**(gamma - 1)*(M - P0)*exp(delta*(-t**gamma + t0**gamma))

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Hyperbolastic_functions
    .. [2] https://tbiomed.biomedcentral.com/articles/10.1186/1742-4682-2-14

    """

    def fdiff(self, argindex=1):
        if argindex == 6:
            M, delta, gama, theta, t0, P0 , t= self.args
            p_t = h3(M, delta, gama, theta, t0, P0 , t)
            return (M - p_t)*(delta*gama*t**(gama-1) + theta/sqrt(1 + \
                                                    theta**2*t**2)).simplify()
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls ,M, delta, gama, theta, t0, P0 , t):
            if sympify(t).is_real and t <= 0 :
                return ValueError("t must be greater than 0")
            else:

                alpha = (M - P0)*exp(delta*t0**gama + asinh(theta*t0))
                return M - alpha*exp(-delta*t**gama-asinh(theta*t)).simplify()
    