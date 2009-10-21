
from sympy.core.basic import Basic, S, C, sympify
from sympy.core.function import Lambda, Function

from sympy.core.evalf import get_integer_part, PrecisionExhausted
from sympy.utilities.decorator import deprecated

###############################################################################
######################### FLOOR and CEILING FUNCTIONS #########################
###############################################################################

class RoundFunction(Function):

    nargs = 1

    @classmethod
    @deprecated
    def canonize(cls, arg):
        return cls.eval(arg)

    @classmethod
    def eval(cls, arg):
        if arg.is_integer:
            return arg
        if arg.is_imaginary:
            return cls(C.im(arg))*S.ImaginaryUnit

        v = cls._eval_number(arg)
        if v is not None:
            return v

        # Integral, numerical, symbolic part
        ipart = npart = spart = S.Zero

        # Extract integral (or complex integral) terms
        if arg.is_Add:
            terms = arg.args
        else:
            terms = [arg]

        for t in terms:
            if t.is_integer or (t.is_imaginary and C.im(t).is_integer):
                ipart += t
            elif t.atoms(C.Symbol):
                spart += t
            else:
                npart += t

        if not (npart or spart):
            return ipart

        # Evaluate npart numerically if independent of spart
        orthogonal = (npart.is_real and spart.is_imaginary) or \
            (npart.is_imaginary and spart.is_real)
        if npart and ((not spart) or orthogonal):
            try:
                re, im = get_integer_part(npart, cls._dir, {}, return_ints=True)
                ipart += C.Integer(re) + C.Integer(im)*S.ImaginaryUnit
                npart = S.Zero
            except (PrecisionExhausted, NotImplementedError):
                pass

        spart = npart + spart
        if not spart:
            return ipart
        elif spart.is_imaginary:
            return ipart + cls(C.im(spart),evaluate=False)*S.ImaginaryUnit
        else:
            return ipart + cls(spart, evaluate=False)

    def _eval_is_bounded(self):
        return self.args[0].is_bounded

    def _eval_is_real(self):
        return self.args[0].is_real

    def _eval_is_integer(self):
        return self.args[0].is_real

class floor(RoundFunction):
    """
    Floor is a univariate function which returns the largest integer
    value not greater than its argument. However this implementaion
    generalizes floor to complex numbers.

    More information can be found in "Concrete mathematics" by Graham,
    pp. 87 or visit http://mathworld.wolfram.com/FloorFunction.html.

        >>> from sympy import *
        >>> floor(17)
        17
        >>> floor(Rational(23, 10))
        2
        >>> floor(2*E)
        5
        >>> floor(-Real(0.567))
        -1
        >>> floor(-I/2)
        -I

    """
    _dir = -1

    @classmethod
    def _eval_number(cls, arg):
        if arg.is_Number:
            if arg.is_Rational:
                if not arg.q:
                    return arg
                return C.Integer(arg.p // arg.q)
            elif arg.is_Real:
                return C.Integer(int(arg.floor()))
        if arg.is_NumberSymbol:
            return arg.approximation_interval(C.Integer)[0]

    def _eval_nseries(self, x, x0, n):
        r = self.subs(x, x0)
        args = self.args[0]
        if args.subs(x, x0) == r:
            direction = (args.subs(x, x+x0) - args.subs(x, x0)).leadterm(x)[0]
            if direction.is_positive:
                return r
            else:
                return r-1
        else:
            return r


class ceiling(RoundFunction):
    """
    Ceiling is a univariate function which returns the smallest integer
    value not less than its argument. Ceiling function is generalized
    in this implementation to complex numbers.

    More information can be found in "Concrete mathematics" by Graham,
    pp. 87 or visit http://mathworld.wolfram.com/CeilingFunction.html.

        >>> from sympy import *
        >>> ceiling(17)
        17
        >>> ceiling(Rational(23, 10))
        3
        >>> ceiling(2*E)
        6
        >>> ceiling(-Real(0.567))
        0
        >>> ceiling(I/2)
        I

    """
    _dir = 1

    @classmethod
    def _eval_number(cls, arg):
        if arg.is_Number:
            if arg.is_Rational:
                if not arg.q:
                    return arg
                return -C.Integer(-arg.p // arg.q)
            elif arg.is_Real:
                return C.Integer(int(arg.ceiling()))
        if arg.is_NumberSymbol:
            return arg.approximation_interval(C.Integer)[1]

    def _eval_nseries(self, x, x0, n):
        r = self.subs(x, x0)
        args = self.args[0]
        if args.subs(x,x0) == r:
            direction = (args.subs(x, x+x0) - args.subs(x, x0)).leadterm(x)[0]
            if direction.is_positive:
                return r+1
            else:
                return r
        else:
            return r

