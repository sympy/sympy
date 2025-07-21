from sympy.core.symbol import symbols
from sympy.calculus.util import Derivative
from sympy.simplify.simplify import simplify
from sympy.integrals.transforms import inverse_laplace_transform
from sympy.simplify.radsimp import apart
from sympy.abc import t, s

def laplace_solve(ode, y, y0=0, y1=0):
    """
    Solve linear ODEs using the Laplace transform method.

    Parameters
    ----------
    ode : Eq
        The differential equation to solve (sympy Eq).
    y : Function
        The dependent function y(t).
    y0 : number, optional
        Initial condition y(0).
    y1 : number, optional
        Initial condition y'(0), used for second-order ODEs.

    Returns
    -------
    Expr
        The time-domain solution y(t).
    """
    y_laplace = symbols('Y')  # Y(s)
    derivatives = {
        Derivative(y, t): s*y_laplace - y0,
        Derivative(y, t, 2): s**2*y_laplace - s*y0 - y1,
    }

    # Bring everything to LHS
    expr = ode.lhs - ode.rhs
    expr = expr.subs(derivatives)
    expr = expr.subs(y, y_laplace)

    Y_s = solve(expr, y_laplace)[0]
    Y_s = apart(simplify(Y_s))

    y_t = inverse_laplace_transform(Y_s, s, t)
    return simplify(y_t)



    



    

