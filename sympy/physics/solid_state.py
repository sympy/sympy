from __future__ import print_function, division

from sympy import S, coth


def brillouin(J, x):
    """
    Returns the Brillouin Function B_{J}(x) that describes the behaviour
    of a paramagnets magnetization in an external magnetic field

    J
        total angular momentum quantum number
    x
        ratio of the energy of the Bohr magneton mu_B inside an external
        magnetic field with flux density B and the thermal energy k_B T;
        multiplied with both the Lande factor g and the total angular
        momentum quantum number J
        x=g*J*B*mu_B/(k_B*T)

    Example
    =======

    >>> from sympy.physics.solid_state import brillouin
    >>> from sympy import var
    >>> var("x J")
    (x, J)
    >>> brillouin(J, x)
    (2*J + 1)*coth(x*(2*J + 1)/(2*J))/(2*J) - coth(x/(2*J))/(2*J)
    """
    # sympify arguments
    J, x= S(J), S(x)
    # two factors inside the Brillioun function
    a=(2*J+1)/(2*J)
    b=1/(2*J)
    #return Brillouin function
    return a*coth(a*x)-b*coth(b*x)

def langevin(x):
    """
    Returns the Langevin Function L(x) that describes the behaviour of
    a paramagnets magnetization in the classical limit J\rightarrow oo.

    Example
    =======

    >>> from sympy.physics.solid_state import langevin
    >>> from sympy import var
    >>> var("x")
    x
    >>> langevin(x)
    coth(x) - 1/x
    """
    # sympify arguments
    x=S(x)
    #return Langevin function
    return coth(x)-S(1)/x
