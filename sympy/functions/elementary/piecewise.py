from sympy.core import Add, Basic, Expr, Function, S, sympify
from sympy.core.containers import Tuple
from sympy.core.relational import Equality, Relational
from sympy.logic.boolalg import Boolean

from warnings import warn
from sympy.core.compatibility import SymPyDeprecationWarning

class Piecewise(Function):
    """
    Represents a piecewise function.

    Usage:

    ``Piecewise( (expr,cond), (expr,cond), ..., (expr, cond), [otherwise])``

        - Each expr/cond pair is a 2-tuple defining a expression and condition
        - The otherwise expression is optional final argument. If no argument
          is given, the default value is NaN.

    Piecewise also accepts an option ``evaluate``, which is False by default.
    When booleans are passed as conditions or other Piecewise functions are
    passed as expressions, evaluation is automatically enabled. On evaluation:

        - Nested Piecewise functions are simplified (currently only simple
          cases can be handled).
        - Explicitly False conditions are removed.
        - Evaluation stops when at an explicitly True condition. The
          corresponding expression becomes the new otherwise expression.
        - If the only remaining argument is the otherwise expression, the
          otherwise expression is returned.

    Examples
    ========

    Create a Piecewise function and evaluate it at 1 and 5:

        >>> from sympy import Piecewise, log
        >>> from sympy.abc import x
        >>> f = x**2
        >>> g = log(x)
        >>> p = Piecewise( (0, x<-1), (f, x<=1), g)
        >>> p.subs(x,1)
        1
        >>> p.subs(x,5)
        log(5)

    To create a cond checking if a variable lies in an Interval, you can use
    the ``.contains()`` or ``.as_relational()`` methods. Evaluate a Piecewise
    function both inside and outside of the defined intervals:

        >>> from sympy import Interval
        >>> i1 = Interval(-2,0)
        >>> i2 = Interval(0,2)
        >>> p = Piecewise( (-x, i1.contains(x)), (x, i2.as_relational(x)))
        >>> p.subs(x,-1)
        1
        >>> p.subs(x,3)
        nan

    See Also
    ========

    piecewise_fold
    """

    nargs = None
    is_Piecewise = True

    def __new__(cls, *args, **options):
        ecs = [ Tuple(*sympify(arg)) for arg in args if hasattr(arg, '__iter__') ]
        oth = [ arg for arg in args if not hasattr(arg, '__iter__') ]
        evaluate = options.pop('evaluate', False)
        if any([ len(ec) != 2 for ec in ecs ]):
            raise ValueError("Piecewise conditions must be (expr, cond) pairs.")

        # deprecated otherwise behavior
        if len(ecs) > 0 and isinstance(ecs[-1][1], bool) and ecs[-1][1] and len(oth) == 0:
            warn('Pass "otherwise" expression as parameter rather than (expr, True). ' \
                'Adding an otherwise parameter, such as the default NaN, silences this warning.', \
                SymPyDeprecationWarning)

        if len(oth) == 0:
            oth = S.NaN
        elif len(oth) == 1:
            oth = sympify(oth[0])
            if list(args).index(oth) != len(args) - 1:
                raise ValueError("Otherwise must be specified as final argument, found in position %d of %d" \
                    % (list(args).index(oth)+1, len(args)))
        else:
            oth = [ str(o) for o in oth ]
            raise ValueError("Only one otherwise statement can be specified, got: %s" % ', '.join(oth))
        new_args = ecs + [oth]

        # evaluation
        if not evaluate and any([ isinstance(expr, Piecewise) or isinstance(cond, bool) for (expr, cond) in ecs ]):
            evaluate = True
        if evaluate:
            evaluated = cls.eval(*new_args)
            if evaluated is not None:
                return evaluated
        # check exprs and otherwise are Exprs
        # TODO: UndefinedFunction neither Basic not Expr
        #from sympy.geometry.entity import GeometryEntity
        #if not all([ isinstance(expr, Expr) or isinstance(expr, GeometryEntity) for (expr, _) in ecs ]):
        #    bad_args = [ "%s of type %s" % (expr, type(expr)) for (expr, _) in ecs if not isinstance(expr, Expr) ]
        #    raise TypeError("Expressions must be subclass of Expr, " \
        #                    "got: %s" % ', '.join(bad_args))
        # check conds are Relational or Boolean
        if not all([ isinstance(cond, Relational) or isinstance(cond, Boolean) for (_, cond) in ecs ]):
            bad_args = [ "%s of type %s" % (cond, type(cond)) for (_, cond) in ecs if not (isinstance(cond, Relational) or isinstance(cond, Boolean)) ]
            raise TypeError("Conditions can only be Relationals or Booleans, " \
                            "got: %s" % ', '.join(bad_args))
        return Expr.__new__(cls, *new_args)

    @classmethod
    def eval(cls, *args):
        from sympy import Or
        # Check for situations where we can evaluate the Piecewise object.
        # 1) Hit an unevaluable cond (e.g. x<1) -> keep object
        # 2) Hit a true condition -> return that expr
        # 3) Remove false conditions, if no conditions left return otherwise
        all_conds_evaled = True    # False if any previous conditions cannot be determined
        piecewise_again = False    # True when args are changed and we should return new Piecewise
        new_args = []
        ecs = args[:-1]
        oth = args[-1]
        or1 = Or( *[c for (_, c) in ecs] )
        for expr, cond in ecs:
            # Check here if expr is a Piecewise and collapse if one of
            # the conds in expr matches cond. This allows the collapsing
            # of Piecewise((Piecewise(x,x<0),x<0)) to Piecewise((x,x<0)).
            # This is important when using piecewise_fold to simplify
            # multiple Piecewise instances having the same conds.
            # Eventually, this code should be able to collapse Piecewise's
            # having different intervals, but this will probably require
            # using the new assumptions.
            if isinstance(expr, Piecewise):
                or2 = Or( *[c for (_, c) in expr.exprcondpairs] )
                for e, c in expr.exprcondpairs:
                    # Don't collapse if cond is "True" as this leads to
                    # incorrect simplifications with nested Piecewises.
                    if c == cond and (or1 == or2 or cond is not True):
                        expr = e
                        piecewise_again = True
            cond_eval = cls.__eval_cond(cond)
            if cond_eval is None:
                all_conds_evaled = False
                new_args.append( (expr, cond) )
            elif cond_eval:
                if all_conds_evaled:
                    return expr
                oth = expr
                break
        if isinstance(oth, Piecewise):
            or2 = Or( *[c for (_, c) in oth.exprcondpairs] )
            if or1 == or2:
                piecewise_again = True
                oth = oth.otherwise
        new_args.append(oth)
        if len(new_args) == 1:
            return oth
        if len(new_args) != len(args) or piecewise_again:
            return Piecewise(*new_args)
        return None

    @property
    def exprcondpairs(self):
        return self.args[:-1]

    @property
    def otherwise(self):
        return self.args[-1]

    @property
    def is_commutative(self):
        exprs = [expr for (expr,_) in self.exprcondpairs]
        oth = self.otherwise
        return all(expr.is_commutative for expr in exprs) and oth.is_commutative

    def _eval_derivative(self, s):
        from sympy.core.function import diff
        new_args = [ (diff(e, s), c) for e, c in self.exprcondpairs ]
        # TODO: diff(S.NaN, s) == 0
        #new_args.append( diff(self.otherwise, s) )
        if self.otherwise is not S.NaN:
            new_args.append( diff(self.otherwise, s) )
        else:
            new_args.append( S.NaN )
        return self.func(*new_args)

    def _eval_integral(self,x):
        from sympy.integrals import integrate
        new_args = [(integrate(e, x), c) for e, c in self.exprcondpairs]
        new_args.append(integrate(self.otherwise, x))
        return self.func(*new_args)

    def _eval_interval(self, sym, a, b):
        """Evaluates the function along the sym in a given interval ab"""
        # FIXME: Currently complex intervals are not supported.  A possible
        # replacement algorithm, discussed in issue 2128, can be found in the
        # following papers;
        #     http://portal.acm.org/citation.cfm?id=281649
        #     http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.70.4127&rep=rep1&type=pdf
        int_expr = []
        mul = 1
        if a > b:
            a, b, mul = b, a, -1
        default = self.otherwise

        # Determine what intervals the expr,cond pairs affect.
        # 1) If cond is True, then log it as default
        # 1.1) Currently if cond can't be evaluated, throw NotImplementedError.
        # 2) For each inequality, if previous cond defines part of the interval
        #    update the new conds interval.
        #    -  eg x < 1, x < 3 -> [oo,1],[1,3] instead of [oo,1],[oo,3]
        # 3) Sort the intervals to make it easier to find correct exprs
        for expr, cond in self.exprcondpairs:
            if cond is True:
                default = expr
                break
            elif cond is False:
                continue
            elif isinstance(cond, Equality):
                continue

            lower, upper = cond.lts, cond.gts # part 1: initialize with givens
            if cond.lts.has(sym):     # part 1a: expand the side ...
                lower = S.NegativeInfinity   # e.g. x <= 0 ---> -oo <= 0
            elif cond.gts.has(sym):   # part 1a: ... that can be expanded
                upper = S.Infinity           # e.g. x >= 0 --->  oo >= 0
            else:
                raise NotImplementedError(
                        "Unable to handle interval evaluation of expression.")

            # part 1b: Reduce (-)infinity to what was passed in.
            lower, upper = max(a, lower), min(b, upper)

            for n in xrange(len(int_expr)):
                # Part 2: remove any interval overlap.  For any conflicts, the
                # interval already there wins, and the incoming interval updates
                # its bounds accordingly.
                if self.__eval_cond(lower < int_expr[n][1]) and \
                        self.__eval_cond(lower >= int_expr[n][0]):
                    lower = int_expr[n][1]
                if self.__eval_cond(upper > int_expr[n][0]) and \
                        self.__eval_cond(upper <= int_expr[n][1]):
                    upper = int_expr[n][0]
            if self.__eval_cond(lower < upper):  # Is it still an interval?
                int_expr.append((lower, upper, expr))
        int_expr.sort(key=lambda x:x[0])

        # Add holes to list of intervals if there is a default value,
        # otherwise raise a ValueError.
        holes = []
        curr_low = a
        for int_a, int_b, expr in int_expr:
            if curr_low < int_a:
                holes.append([curr_low, min(b, int_a), default])
            curr_low = int_b
            if curr_low > b:
                break
        if curr_low < b:
            holes.append([curr_low, b, default])

        if holes:
            int_expr.extend(holes)

        # Finally run through the intervals and sum the evaluation.
        ret_fun = 0
        for int_a, int_b, expr in int_expr:
            ret_fun += expr._eval_interval(sym,  max(a, int_a), min(b, int_b))
        return mul * ret_fun

    def _eval_nseries(self, x, n, logx):
        new_args = map(lambda ec: (ec[0]._eval_nseries(x, n, logx), ec[1]), self.exprcondpairs)
        new_args.append( self.otherwise._eval_nseries(x, n, logx) )
        return self.func( *new_args )

    def _eval_subs(self, old, new):
        new_args = [(expr._eval_subs(old, new), cond._eval_subs(old, new)) for expr, cond in self.exprcondpairs]
        new_args.append(self.otherwise._eval_subs(old, new))
        return self.func(*new_args, **{'evaluate': True})

    def _eval_as_leading_term(self, x):
        args = [ e for (e, _) in self.exprcondpairs]
        if self.otherwise is not S.NaN:
            args.append( self.otherwise )
        return Add(*args).as_leading_term(x)

    @classmethod
    def __eval_cond(cls, cond):
        """Returns cond if it's a number, otherwise it is undecidable and returns None."""
        if isinstance(cond,bool):
            return cond
        return None

def piecewise_fold(expr):
    """
    Takes an expression containing a piecewise function and returns the
    expression in piecewise form.

    Examples
    ========

    >>> from sympy import Piecewise, piecewise_fold, sympify as S
    >>> from sympy.abc import x
    >>> p = Piecewise((x, x < 1), (1, S(1) <= x))
    >>> piecewise_fold(x*p)
    Piecewise((x**2, x < 1), (x, 1 <= x))

    See Also
    ========

    Piecewise
    """
    if not isinstance(expr, Basic) or not expr.has(Piecewise):
        return expr
    new_args = map(piecewise_fold, expr.args)
    if expr.func is Tuple:
        return Tuple(*new_args)
    piecewise_args = []
    for n, arg in enumerate(new_args):
        if sympify(arg).func is Piecewise:
            piecewise_args.append(n)
    if len(piecewise_args) > 0:
        n = piecewise_args[0]
        ecs = new_args[n].exprcondpairs
        oth = new_args[n].otherwise
        prev_args = new_args[:n]
        post_args = new_args[n+1:]
        new_args = [(expr.func(*(prev_args + [e] + post_args)), c) \
                for e, c in ecs]
        new_args.append(expr.func(*(prev_args + [oth] + post_args)))
        if len(piecewise_args) > 1:
            return piecewise_fold(Piecewise(*new_args))
    return Piecewise(*new_args)
