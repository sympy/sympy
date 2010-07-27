from sympy import factorial, sqrt, exp, S, laguerre_l

def R_nl(n, l, a, r):
    """
    Returns the Hydrogen radial wavefunction R_{nl}.

    n, l .... quantum numbers 'n' and 'l'
    a    .... Bohr radius \hbar^2 / (Z*M)
    r    .... radial coordinate

    Examples::

    >>> from sympy.physics.hydrogen import R_nl
    >>> from sympy import var
    >>> var("r a")
    (r, a)
    >>> R_nl(1, 0, a, r)
    2*(a**(-3))**(1/2)*exp(-r/a)
    >>> R_nl(2, 0, a, r)
    2**(1/2)*(a**(-3))**(1/2)*(2 - r/a)*exp(-r/(2*a))/4
    >>> R_nl(2, 1, a, r)
    r*6**(1/2)*(a**(-3))**(1/2)*exp(-r/(2*a))/(12*a)

    In most cases you probably want to set a=1::

    >>> R_nl(1, 0, 1, r)
    2*exp(-r)
    >>> R_nl(2, 0, 1, r)
    2**(1/2)*(2 - r)*exp(-r/2)/4
    >>> R_nl(3, 0, 1, r)
    2*3**(1/2)*(3 - 2*r + 2*r**2/9)*exp(-r/3)/27

    The normalization of the radial wavefunction is::

    >>> from sympy import integrate, oo
    >>> integrate(R_nl(1, 0, 1, r)**2 * r**2, (r, 0, oo))
    1
    >>> integrate(R_nl(2, 0, 1, r)**2 * r**2, (r, 0, oo))
    1
    >>> integrate(R_nl(2, 1, 1, r)**2 * r**2, (r, 0, oo))
    1

    """
    # sympify arguments
    n, l, a, r = S(n), S(l), S(a), S(r)
    # radial quantum number
    n_r = n - l - 1
    # rescaled "r"
    r0 = 2 * r / (n * a)
    # normalization coefficient
    C =  sqrt((S(2)/(n*a))**3 * factorial(n_r) / (2*n*factorial(n+l)))
    # This is an equivalent normalization coefficient, that can be found in
    # some books. Both coefficients seem to be the same fast:
    # C =  S(2)/n**2 * sqrt(1/a**3 * factorial(n_r) / (factorial(n+l)))
    return  C * r0**l * laguerre_l(n_r, 2*l+1, r0) * exp(-r0/2)
