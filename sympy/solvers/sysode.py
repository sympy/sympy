"""
This module implements solving systems of differential equations
"""

from sympy.core import Add, Basic, C, S, Mul
from sympy.core.relational import Equality, Eq
from sympy.core.symbol import Symbol
from sympy.core.sympify import sympify

from sympy.functions import exp
from sympy.matrices import zeros, Matrix
from sympy.simplify import collect, simplify

from sympy.abc import x,y,z

def ode_1st_linear_system(eqs, wrt, *symbols):
    """
    Solves a system of first order homogenous linear differential equations

    An example of such a system is a set equations such as
    dx/dt = a1*x + b1*y + c1(t), dy/dt = a2*x + b2*y + c2(t) where a1, b1, a2 and b2
    are typically arbitrary scalars and c1 and c2 are functions of t.
    To solve such equations we generally need to find the eigenvalues and
    eigenvectors of the matrix of scalars and construct the functions
    we need from there.

    ***Examples***
        This is the example given in wikipedia
        >>> from sympy.abc import x,y,z
        >>> from sympy.solvers.sde import *
        >>> t = Symbol('t')
        >>> ode_1st_linear_system([3*y -4*z, 4*y - 7*z], t, y, z)
        [y == 2*C0*exp(t) + C1*exp(-5*t)/2, z == C0*exp(t) + C1*exp(-5*t)]

        We can put in arbitrary functions as well
        >>> ode_1st_linear_system([y - z + t, -y + z + 1 ], t, y, z)
        [y == C1 + t/2 + t**2/4 - (C0 - exp(-2*t)/8 + t*exp(-2*t)/4)*exp(2*t),\
        z == C1 + t/2 + t**2/4 + (C0 - exp(-2*t)/8 + t*exp(-2*t)/4)*exp(2*t)]

        
    """
    if not len(symbols) == len(eqs):
        raise ValueError("The number of symbols must equal the number of equations")

    scalar_mat = zeros(len(eqs))
    arbitrary_funcs = zeros((len(eqs), 1))

    symbol_map = {}
    for i in range(len(symbols)):
        symbol_map[symbols[i]] = i

    for i, eq in enumerate(eqs):
        # Iterate over all the equations
        eq = collect(eq, *symbols)
        if not eq.is_Add:
            term_tup = eq.as_independent(*symbols)
            base_term = term_tup[1]
            if base_term in symbols:
                scalar_mat[i, symbol_map[base_term]] = term_tup[0]
            else:
                arbitrary_funcs[i] += eq
        else:
            eq = collect(eq, *symbols)
            for j in eq.args:
                # Iterate over all the terms
                term_tup = j.as_independent(*symbols)
                base_term = term_tup[1]
                if term_tup[0] == j:
                    arbitrary_funcs[i] += j
                    continue

                if base_term in symbols:
                    scalar_mat[i, symbol_map[base_term]] = term_tup[0]
                    
    # We need the eigenvalues and eigenvectors
    eigenvals = scalar_mat.berkowitz_eigenvals()
    eigenval_list = []
    for k, v in eigenvals.iteritems():
        for i in range(v):
            eigenval_list.append(k)

    eigenvects = scalar_mat.eigenvects()
    eigen_mat = zeros((len(eqs), 1))
    for entry in eigenvects:
        entry_eigenvals = entry[2]
        for vect in entry_eigenvals:
            eigen_mat = eigen_mat.row_join(vect)

    eigen_mat.col_del(0)
    arbitrary_funcs = eigen_mat.inv().multiply(arbitrary_funcs)
    consts = [Symbol("C"+str(i)) for i in range(len(eqs))]

    # Construct the function
    funcs = Matrix([S((consts[i] + C.Integral(exp(-eigenval_list[i]*wrt)*arbitrary_funcs[i],wrt).doit())\
                      *exp(eigenval_list[i]*wrt)) for i in range(len(eqs))])
    return [Eq(symbols[i], eigen_mat.multiply(funcs).tolist()[i][0]) for i in range(len(eqs))]
