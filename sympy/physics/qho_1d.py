from sympy.core import S, pi, Rational
from sympy.functions import hermite, sqrt, exp, factorial

def X_n(n, nu, x):
    """
    Returns the wavefunction X_{n} for the One-dimensional harmonic oscillator.

    ``n``
        the "nodal" quantum number.  Corresponds to the number of nodes in the
        wavefunction.  n >= 0
    ``nu``
        mass-scaled frequency: nu = m*omega/(hbar) where `m' is the mass and
        `omega' the frequency of the oscillator.  (in atomic units nu == omega/2)
    ``x``
        x coordinate
        
    :Example:

    >>> from sympy.physics.qho_1d import X_n
    >>> from sympy import var
    >>> var("x nu")
    (x, nu)
    >>> X_n(0, nu, x)
    nu**(1/4)*exp(-nu*x**2/2)/pi**(1/4)
    
    """
    
    # sympify arguments
    n, nu, x = map(S, [n, nu, x])
    # normalization coefficient
    C =  (nu/pi)**(S(1)/4) * sqrt(1/(2**n*factorial(n)))

    return  C * exp(-nu* x**2 /2) * hermite(n, sqrt(nu)*x)
    
    
def E_n(n,hw):
    """
    Returns the Energy of the One-dimensional harmonic oscillator

    ``n``
        the "nodal" quantum number
    ``hw``
        the harmonic oscillator parameter.

    The unit of the returned value matches the unit of hw, since the energy is
    calculated as:

        E_nl = hw*(n + 1/2)
    """
  
    return hw*(n + Rational(1,2))
    