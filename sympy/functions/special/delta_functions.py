from sympy.core import S, sympify, Function, diff
import sympy.polys

###############################################################################
################################ DELTA FUNCTION ###############################
###############################################################################
class DiracDelta(Function):
    """DiracDelta function, and the derivatives.
    DiracDelta function has the following properties:
    1) diff(Heaviside(x),x) = DiracDelta(x)
    2) integrate(DiracDelta(x-a)*f(x),(x,-oo,oo)) = f(a)
       integrate(DiracDelta(x-a)*f(x),(x,a-e,a+e)) = f(a)
    3) DiracDelta(x) = 0, for all x != 0
    4) DiracDelta(g(x)) = Sum_i(DiracDelta(x-xi)/abs(g'(xi)))
       Where xis are the roots of g

    Derivatives of k order of DiracDelta have the following property:
    5) DiracDelta(x,k) = 0, for all x!=0


    For more information, see:
    http://mathworld.wolfram.com/DeltaFunction.html
    """

    nargs = (1,2)

    def fdiff(self, argindex = 1):
        if argindex == 1:
            #I didn't know if there is a better way to handle default arguments
            k = 0
            if len(self.args) > 1:
                k = self.args[1]
            return DiracDelta(self.args[0],k+1)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg, k=0):
        k = sympify(k)
        if not k.is_Integer or k.is_negative:
            raise ValueError("Error: the second argument of DiracDelta must be \
            a non-negative integer, %s given instead." %(k,))
        arg = sympify(arg)
        if arg is S.NaN:
            return S.NaN
        if arg.is_positive or arg.is_negative:
            return S.Zero
        elif arg.is_zero:
            return S.Infinity

    def simplify(self, x):
        """simplify(self, x)

           Compute a simplified representation of the function using
           property number 4.

           x can be:

           - a symbol

           Examples
           --------

           >>> from sympy import DiracDelta
           >>> from sympy.abc import x, y

           >>> DiracDelta(x*y).simplify(x)
           DiracDelta(x)/abs(y)
           >>> DiracDelta(x*y).simplify(y)
           DiracDelta(y)/abs(x)

           >>> DiracDelta(x**2+x-2).simplify(x)
           DiracDelta(-1 + x)/3 + DiracDelta(2 + x)/3

        """
        if not self.args[0].has(x) or (len(self.args)>1 and self.args[1] != 0 ):
            return self
        try:
            argroots = sympy.polys.polyroots.roots(self.args[0],x, \
                                                     multiple=True)
            result = 0
            valid = True
            darg = diff(self.args[0], x)
            for r in argroots:
                #should I care about multiplicities of roots?
                if r.is_real and not darg.subs(x,r).is_zero:
                    result = result + DiracDelta(x - r)/abs(darg.subs(x,r))
                else:
                    valid = False
                    break
            if valid:
                return result
        except Exception,e :
            print e
            raise
            pass
        return self

    def is_simple(self,x):
        """is_simple(self, x)

           Tells whether the argument(args[0]) of DiracDelta is a linear
           expression in x.

           x can be:

           - a symbol

           Examples
           --------

           >>> from sympy import DiracDelta, cos
           >>> from sympy.abc import x, y

           >>> DiracDelta(x*y).is_simple(x)
           True
           >>> DiracDelta(x*y).is_simple(y)
           True

           >>> DiracDelta(x**2+x-2).is_simple(x)
           False

           >>> DiracDelta(cos(x)).is_simple(x)
           False

        """
        p = self.args[0].as_poly(x)
        if p:
            return p.degree() == 1
        return False
###############################################################################
############################## HEAVISIDE FUNCTION #############################
###############################################################################

class Heaviside(Function):
    """Heaviside Piecewise function.
    Heaviside function has the following properties:
    1) diff(Heaviside(x),x) = DiracDelta(x)
                        ( 0, if x<0
    2) Heaviside(x) = < [*]  1/2 if x==0
                        ( 1, if x>0
    [*]Regarding to the value at 0, Mathematica adopt the value H(0)=1,
    and Maple H(0)=undefined

    I think is better to have H(0)=1/2, due to the following:
    integrate(DiracDelta(x),x) = Heaviside(x)
    integrate(DiracDelta(x),(x,-oo,oo)) = 1

    and since DiracDelta is a symmetric function,
    integrate(DiracDelta(x),(x,0,oo)) should be 1/2
    in fact, that is what maple returns.

    If we take Heaviside(0)=1/2, we would have
    integrate(DiracDelta(x),(x,0,oo)) = Heaviside(oo)-Heaviside(0)=1-1/2= 1/2
    and
    integrate(DiracDelta(x),(x,-oo,0)) = Heaviside(0)-Heaviside(-oo)=1/2-0= 1/2

    If we consider, instead Heaviside(0)=1, we would have
    integrate(DiracDelta(x),(x,0,oo)) = Heaviside(oo)-Heaviside(0) = 0
    and
    integrate(DiracDelta(x),(x,-oo,0)) = Heaviside(0)-Heaviside(-oo) = 1


    For more information, see:
    http://mathworld.wolfram.com/HeavisideStepFunction.html
    """
    nargs = 1

    def fdiff(self, argindex = 1):
        if argindex == 1:
            # property number 1
            return DiracDelta(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        arg = sympify(arg)
        if arg is S.NaN:
            return S.NaN
        elif arg.is_negative:
            return S.Zero
        elif arg.is_zero:
            return S.Half
        elif arg.is_positive:
            return S.One

