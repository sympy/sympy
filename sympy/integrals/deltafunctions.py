import sympy
from sympy.core import Symbol, Wild, S
from sympy.functions import DiracDelta, Heaviside
from sympy.solvers import solve
#from sympy.integrals import Integral

def change_mul(node,x):
    """change_mul(node,x)

       Rearranges the operands of a product, bringing to front any simple
       DiracDelta expression.

       If no simple DiracDelta expression was found, then all the DiracDelta
       expressions are simplified (using DiracDelta.simplify).

       Return: (dirac,nnode)
       Where:
       dirac is a simple DiracDelta expression. None if no simple expression has been found
       nnode is a new node where all the DiracDelta expressions where simplified,
       and finally the node was expanded. if nnode is None, means that no DiracDelta expression
       could be simplified

       Examples
       --------

       >>change_mul(x*y*DiracDelta(x)*cos(x),x)
       (DiracDelta(x),x*y*cos(x))
       >>change_mul(x*y*DiracDelta(x**2-1)*cos(x),x)
       (None,x*y*cos(x),x*y*DiracDelta(1 + x)*cos(x)/2 + x*y*DiracDelta(-1 + x)*cos(x)/2)
       >>change_mul(x*y*DiracDelta(cos(x))*cos(x),x)
       (None,None)

    """
    if not node.is_Mul:
        return node
    new_args = []
    dirac = None
    for arg in node.args:
        if arg.func == DiracDelta and arg.is_simple(x) \
                and (len(arg.args) <= 1 or arg.args[1]==0):
            dirac = arg
        else:
            new_args.append(change_mul(arg,x))
    if not dirac:#we didn't find any simple dirac
        new_args = []
        for arg in node.args:
            if arg.func == DiracDelta:
                new_args.append(arg.simplify(x))
            else:
                new_args.append(change_mul(arg,x))
        if tuple(new_args) != node.args:
            nnode = node.__class__(*new_args).expand()
        else:#if the node didn't change there is nothing to do
            nnode = None
        return (None, nnode)
    return (dirac, node.func(*new_args))


def deltaintegrate(f, x):
    """The idea for integration is the following:
    -If we are dealing with a DiracDelta expression, i.e.:
    DiracDelta(g(x)), we try to simplify it.
    If we could simplify it, then we integrate the resulting expression.
    We already know we can integrate a simplified expression, because only
    simple DiracDelta expressions are involved.
    If we couldn't simplify it, there are two cases:
    1) The expression is a simple expression, then we return the integral
    Taking care if we are dealing with a Derivative or with a proper DiracDelta
    2) The expression is not simple(i.e. DiracDelta(cos(x))), we can do nothing at all

    -If the node is a multiplication node having a DiracDelta term
    First we expand it.
    If the expansion did work, the we try to integrate the expansion
    If not, we try to extract a simple DiracDelta term, then we have two cases
    1)We have a simple DiracDelta term, so we return the integral
    2)We didn't have a simple term, but we do have an expression with simplified
    DiracDelta terms, so we integrate this expression

    """
    if not f.has(DiracDelta):
        return None
    # g(x) = DiracDelta(h(x))
    if f.func == DiracDelta:
        h = f.simplify(x)
        if h == f:#can't simplify the expression
            #FIXME: the second term tells whether is DeltaDirac or Derivative
            #For integrating derivatives of DiracDelta we need the chain rule
            if f.is_simple(x):
                if (len(f.args) <= 1 or f.args[1]==0):
                    return Heaviside(f.args[0])
                else:
                    return (DiracDelta(f.args[0],f.args[1]-1)/ f.args[0].as_poly().coeffs[0])
        else:#let's try to integrate the simplified expression
            fh = sympy.integrals.integrate(h,x)
            return fh
    elif f.is_Mul: #g(x)=a*b*c*f(DiracDelta(h(x)))*d*e
        g = f.expand()
        if f != g:#the expansion worked
            fh = sympy.integrals.integrate(g,x)
            if fh and not isinstance(fh,sympy.integrals.Integral):
                return fh
        else:#no expansion performed, try to extract a simple DiracDelta term
            dg, rest_mult = change_mul(f,x)
            if not dg:
                if rest_mult:
                    fh = sympy.integrals.integrate(rest_mult,x)
                    return fh
            else:
                point = solve(dg.args[0],x)[0]
                return (rest_mult.subs(x,point)*Heaviside(dg.args[0]))
    return None
