
from sympy.core.basic import Basic, S
from sympy.core.cache import cacheit
from sympy.core.function import Function, FunctionClass, diff
from sympy.core.numbers import Number
from sympy.core.relational import Relational
from sympy.core.sympify import sympify

class Piecewise(Function):
    """
    Represents a piecewise function.

    Usage
    =====
      Piecewise( (expr,cond), (expr,cond), ... )
        - Each argument is a 2-tuple defining a expression and condition
        - The conds are evaluated in turn returning the first that is True.
          If any of the evaluated conds are not determined explicitly False,
          e.g. x < 1, the function is returned in symbolic form.
        - If the function is evaluated at a place where all conditions are False,
          a ValueError exception will be raised.
        - Pairs where the cond is explicitly False, will be removed.

    Examples
    ========
      >>> from sympy import *
      >>> x = Symbol('x')
      >>> f = x**2
      >>> g = log(x)
      >>> p = Piecewise( (0, x<-1), (f, x<=1), (g, True))
      >>> p.subs(x,1)
      1
      >>> p.subs(x,5)
      log(5)
    """

    nargs=None

    @cacheit
    def __new__(cls, *args, **options):
        for opt in ["nargs", "dummy", "comparable", "noncommutative", "commutative"]:
            if opt in options:
                del options[opt]
        r = cls.canonize(*args)
        if r is None:
            return Basic.__new__(cls, *args, **options)
        else:
            return r

    @classmethod
    def canonize(cls, *args):
        # Check types first
        for ec in args:
            if (not isinstance(ec,tuple)) or len(ec)!=2:
                raise TypeError, "args may only include (expr, cond) pairs"
        for expr, cond in args:
            cond_type = type(ec[1])
            if cond_type != bool and not issubclass(cond_type,Relational):
                raise TypeError, \
                    "Cond %s is of type %s, but must be a bool or Relational" \
                    % (cond, cond_type)

        # Check for situations where we can evaluate the Piecewise object.
        # 1) Hit an unevaluatable cond (e.g. x<1) -> keep object
        # 2) Hit a true condition -> return that expr
        # 3) Remove false conditions, if no conditions left -> raise ValueError
        all_conds_evaled = True
        non_false_ecpairs = []
        for expr, cond in args:
            cond_eval = cls.__eval_cond(cond)
            if cond_eval is None:
                all_conds_evaled = False
                non_false_ecpairs.append( (expr, cond) )
            elif cond_eval:
                if all_conds_evaled:
                    return expr
                non_false_ecpairs.append( (expr, cond) )
        if len(non_false_ecpairs) != len(args):
            return Piecewise(*non_false_ecpairs)

        # Count number of arguments.
        nargs = 0
        for expr, cond in args:
            if hasattr(expr, 'nargs'):
                nargs = max(nargs, expr.nargs)
            elif hasattr(expr, 'args'):
                nargs = max(nargs, len(expr.args))
        if nargs:
            cls.narg = nargs
        return None

    def doit(self, **hints):
        new_ecpairs = []
        for expr, cond in self.args:
            if hasattr(expr,'doit'):
                new_expr = expr.doit(**hints)
            else:
                new_expr = expr
            if hasattr(cond,'doit'):
                new_cond = cond.doit(**hints)
            else:
                new_cond = cond
            new_ecpairs.append( (new_expr, new_cond) )
        return Piecewise(*new_ecpairs)

    def _eval_integral(self,x):
        from sympy.integrals import integrate
        return  Piecewise(*[(integrate(expr, x),cond) \
                                for expr, cond in self.args])

    def _eval_interval(self, sym, ab):
        """Evaluates the function along the sym in a given interval ab"""
        # FIXME: Currently only supports conds of type sym < Num, or Num < sym
        int_expr = []
        a, b = ab
        mul = 1
        if a > b:
            a = ab[1]; b = ab[0]; mul = -1
        default = None

        # Determine what intervals the expr,cond pairs affect.
        # 1) If cond is True, then log it as default
        # 1.1) Currently if cond can't be evaluated, throw NotImplentedError.
        # 2) For each inequality, if previous cond defines part of the interval
        #    update the new conds interval.
        #    -  eg x < 1, x < 3 -> [oo,1],[1,3] instead of [oo,1],[oo,3]
        # 3) Sort the intervals to make it easier to find correct exprs
        for expr, cond in self.args:
            if type(cond) == bool:
                if cond:
                    default = expr
                    break
                else:
                    continue
            curr = list(cond.args)
            if cond.args[0] == sym:
                curr[0] = S.NegativeInfinity
            elif cond.args[1] == sym:
                curr[1] = S.Infinity
            else:
                raise NotImplementedError, \
                    "Currently only supporting evaluation with only "\
                    "sym on one side fo the relation."
            curr = [max(a,curr[0]),min(b,curr[1])]
            for n in xrange(len(int_expr)):
                if self.__eval_cond(curr[0] < int_expr[n][1]) and \
                        self.__eval_cond(curr[0] >= int_expr[n][0]):
                    curr[0] = int_expr[n][1]
                if self.__eval_cond(curr[1] > int_expr[n][0]) and \
                        self.__eval_cond(curr[1] <= int_expr[n][1]):
                    curr[1] = int_expr[n][0]
            if self.__eval_cond(curr[0] < curr[1]):
                int_expr.append(curr+[expr])
        int_expr.sort(lambda x,y:1 if x[0] > y[0] else -1)

        # Add holes to list of intervals if there is a default value,
        # otherwise raise a ValueError.
        holes = []
        curr_low = a
        for int_a, int_b, expr in int_expr:
            if curr_low < int_a:
                holes.append([curr_low, min(b,int_a), default])
            curr_low = int_b
            if curr_low > b:
                break
        if holes and default != None:
            int_expr.extend(holes)
        elif holes and default == None:
            raise ValueError, "Called interval evaluation over piecewise "\
                              "function on undefined intervals %s" %\
                              ", ".join([str((h[0],h[1])) for h in holes])

        # Finally run through the intervals and sum the evaluation.
        ret_fun = 0
        for int_a, int_b, expr in int_expr:
            B = expr.subs(sym, min(b,int_b))
            if B is S.NaN:
                B = limit(expr, sym, min(b,int_b))
            if B is S.NaN:
                return self
            A = expr.subs(sym, max(a,int_a))
            if A is S.NaN:
                A = limit(expr, sym, max(a,int_a))
            if A is S.NaN:
                return self
            ret_fun += B - A
        return mul * ret_fun

    def _eval_derivative(self, s):
        return Piecewise(*[(diff(expr,s),cond) for expr, cond in self.args])

    def _eval_subs(self, old, new):
        if self == old:
            return new
        new_ecpairs = []
        for expr, cond in self.args:
            if hasattr(expr,"subs") and not isinstance(expr,FunctionClass):
                new_expr = expr.subs(old, new)
            else:
                new_expr = expr
            if hasattr(cond,"subs"):
                new_cond = cond.subs(old, new)
            else:
                new_cond = cond
            new_ecpairs.append( (new_expr, new_cond) )
        return Piecewise(*new_ecpairs)

    @classmethod
    def __eval_cond(cls, cond):
        """Returns if the condition is True or False.
           If it is undeterminable, returns None."""
        if type(cond) == bool:
            return cond
        arg0 = cond.args[0]
        arg1 = cond.args[1]
        if isinstance(arg0, FunctionClass) or isinstance(arg1, FunctionClass):
            return None
        if hasattr(arg0,'evalf'):
            arg0 = arg0.evalf()
        if not issubclass(type(arg0),Number) and \
                type(arg0) != int and type(arg0) != float:
            return None
        if hasattr(arg1,'evalf'):
            arg1 = arg1.evalf()
        if not issubclass(type(arg1),Number) and \
                type(arg1) != int and type(arg1) != float:
            return None
        return bool(cond)
