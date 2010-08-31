from sympy import factorial, sqrt, exp, S, laguerre_l

def R_nl(n, l, r, Z=1):
    """
    Returns the Hydrogen radial wavefunction R_{nl}.

    n, l .... quantum numbers 'n' and 'l'
    r    .... radial coordinate
    Z    .... atomic number (1 for Hydrogen, 2 for Helium, ...)

    Everything is in Hartree atomic units.

    Examples::

    >>> from sympy.physics.hydrogen import R_nl
    >>> from sympy import var
    >>> var("r Z")
    (r, Z)
    >>> R_nl(1, 0, r, Z)
    2*(Z**3)**(1/2)*exp(-Z*r)
    >>> R_nl(2, 0, r, Z)
    2**(1/2)*(Z**3)**(1/2)*(2 - Z*r)*exp(-Z*r/2)/4
    >>> R_nl(2, 1, r, Z)
    Z*r*6**(1/2)*(Z**3)**(1/2)*exp(-Z*r/2)/12

    For Hydrogen atom, you can just use the default value of Z=1::

    >>> R_nl(1, 0, r)
    2*exp(-r)
    >>> R_nl(2, 0, r)
    2**(1/2)*(2 - r)*exp(-r/2)/4
    >>> R_nl(3, 0, r)
    2*3**(1/2)*(3 - 2*r + 2*r**2/9)*exp(-r/3)/27

    For Silver atom, you would use Z=47::

    >>> R_nl(1, 0, r, Z=47)
    94*47**(1/2)*exp(-47*r)
    >>> R_nl(2, 0, r, Z=47)
    47*94**(1/2)*(2 - 47*r)*exp(-47*r/2)/4
    >>> R_nl(3, 0, r, Z=47)
    94*141**(1/2)*(3 - 94*r + 4418*r**2/9)*exp(-47*r/3)/27

    The normalization of the radial wavefunction is::

    >>> from sympy import integrate, oo
    >>> integrate(R_nl(1, 0, r)**2 * r**2, (r, 0, oo))
    1
    >>> integrate(R_nl(2, 0, r)**2 * r**2, (r, 0, oo))
    1
    >>> integrate(R_nl(2, 1, r)**2 * r**2, (r, 0, oo))
    1

    It holds for any atomic number:

    >>> integrate(R_nl(1, 0, r, Z=2)**2 * r**2, (r, 0, oo))
    1
    >>> integrate(R_nl(2, 0, r, Z=3)**2 * r**2, (r, 0, oo))
    1
    >>> integrate(R_nl(2, 1, r, Z=4)**2 * r**2, (r, 0, oo))
    1

    """
    # sympify arguments
    n, l, r, Z = S(n), S(l), S(r), S(Z)
    # radial quantum number
    n_r = n - l - 1
    # rescaled "r"
    a = 1/Z # Bohr radius
    r0 = 2 * r / (n * a)
    # normalization coefficient
    C =  sqrt((S(2)/(n*a))**3 * factorial(n_r) / (2*n*factorial(n+l)))
    # This is an equivalent normalization coefficient, that can be found in
    # some books. Both coefficients seem to be the same fast:
    # C =  S(2)/n**2 * sqrt(1/a**3 * factorial(n_r) / (factorial(n+l)))
    return  C * r0**l * laguerre_l(n_r, 2*l+1, r0) * exp(-r0/2)
