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

    Note that old otherwise syntax of using including an (otherwise, True)
    tuple is currently deprecated. If your final expr/cond is (expr, True),
    you must specify an otherwise expression to avoid having this expr changed
    to an otherwise condition. The defualt value, NaN, is a suitable option.

    Piecewise can also evaluate the function, which is set by the ``evaluate``
    keyword, that is disabled by default. Having explicitly True or False
    conditions will enable evaluation unless the ``evaluate`` option is
    explicitly False. When evaluate is True, it returns the first expression
    with its condition explicitly True.

    When booleans are given as conditions and ``evaluate=False``:

        - Expr/cond pairs where cond is explicitly False are removed.
        - Evaluation stops at an explicitly True condition. The corresponding
          expression becomes the new otherwise expression.
        - If the only remaining argument is the otherwise expression, the
          otherwise expression is returned.

    When Piecewise functions are nested, they will automatically be siplified.
    This behavior can be overridden with the keyword evaluate=False.

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
        evaluate = options.pop('evaluate', None)
        if any([ len(ec) != 2 for ec in ecs ]):
            raise ValueError("Piecewise conditions must be (expr, cond) pairs.")

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

        # deprecated otherwise behavior
        if len(ecs) > 0 and ecs[-1][1] is True and ecs[-1] == args[-1]:
            oth = ecs[-1][0]
            ecs = ecs[:-1]
            warn('Pass "otherwise" expression as parameter rather than (expr, True). ' \
                'Converting final expr/cond to otherwise to maintain old behavior. ' \
                'Manually adding an otherwise parameter, such as NaN, silences this warning and ensures evaluation.', \
                SymPyDeprecationWarning)

        # Collapse nested piecewise expressions
        # Allow explicit evaluate=False to stop evaluation of nested Piecewise
        if evaluate is not False and (any([ isinstance(expr, Piecewise) for (expr, _) in ecs ]) or isinstance(oth, Piecewise)):
            ecs, oth = cls._collapse_piecewise_args(ecs, oth)
        # TODO: See Issue 3025
        # When evaluate=False, remove False conds and move True cond to otherwise
        # Should subs True for Basic True, False for Basic False when evaluate!=False
        if any([ isinstance(cond, bool) for (_, cond) in ecs ]):
            new_ecs = []
            if evaluate is False:
                for e, c in ecs:
                    cond_eval = cls.__eval_cond(c)
                    if cond_eval is False:
                        continue
                    elif cond_eval is True:
                        oth = e
                        break
                    else:
                        new_ecs.append(Tuple(e,c))
                ecs = new_ecs
            else:
                evaluate = True
        # Having no expr/cond pairs should also trigger evaluation unless explicitly False
        if len(ecs) == 0 and evaluate is not False:
            return oth

        new_args = ecs + [oth]
        if evaluate:
            return cls.eval(*new_args)

        # check exprs and otherwise are Exprs
        from sympy.geometry.entity import GeometryEntity
        if not all([ isinstance(expr, Expr) or isinstance(expr, GeometryEntity) for (expr, _) in ecs ]):
            bad_args = [ "%s of type %s" % (expr, type(expr)) for (expr, _) in ecs \
                         if not isinstance(expr, Expr) ]
            raise TypeError("Expressions must be subclass of Expr, " \
                            "got: %s" % ', '.join(bad_args))
        if not (isinstance(oth, Expr) or isinstance(oth, GeometryEntity)):
            raise TypeError("Otherwise expression must be a subclass of Expr, " \
                            "got %s of type %s" % (oth, type(oth)) )
        # check conds are Relational or Boolean
        if not all([ isinstance(cond, Relational) or isinstance(cond, Boolean) for (_, cond) in ecs ]):
            bad_args = [ "%s of type %s" % (cond, type(cond)) for (_, cond) in ecs \
                         if not (isinstance(cond, Relational) or isinstance(cond, Boolean)) ]
            raise TypeError("Conditions can only be Relationals or Booleans, " \
                            "got: %s" % ', '.join(bad_args))

        return Expr.__new__(cls, *new_args)

    def doit(self, **hints):
        evaluate = hints.get('piecewise', True)
        if hints.get('deep', True):
            args = [ arg.doit(**hints) for arg in self.args ]
        else:
            args = self.args
        return self.func(*args, **{'evaluate': evaluate})

    @classmethod
    def eval(cls, *args):
        """
        Evaluate Piecewise function

        Return expression of the first explicitly True condition. If no
        conditions are explicitly True, return otherwise expression.
        """
        ecs = args[:-1]
        oth = args[-1]
        for expr, cond in ecs:
            cond_eval = cls.__eval_cond(cond)
            if cond_eval:
                return expr
        return oth

    @property
    def exprcondpairs(self):
        """
        Return all expressions and conditions as 2-Tuples

        Returns a tuple of 2-Tuples, each formatted as (expr, cond) for each of
        the given expressions and conditions.

            >>> from sympy import Piecewise
            >>> from sympy.abc import x
            >>> p = Piecewise((x, x > 1), (x**2, x > 0), -x)
            >>> p.exprcondpairs
            ((x, x > 1), (x**2, x > 0))

        See Also
        ========

        otherwise
        """
        return self.args[:-1]

    @property
    def otherwise(self):
        """
        Returns otherwise expression

        Returns the expression given when all conditions are False. If no
        expression was specified, this will return NaN.

            >>> from sympy import Piecewise
            >>> from sympy.abc import x
            >>> p = Piecewise((x, x > 1), (x**2, x > 0), -x)
            >>> p.otherwise
            -x
            >>> p = Piecewise((x, x > 1), (x**2, x > 0))
            >>> p.otherwise
            nan

        See Also
        ========

        exprcondpairs
        """
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

    def _eval_as_leading_term(self, x):
        args = [ e for (e, _) in self.exprcondpairs]
        if self.otherwise is not S.NaN:
            args.append( self.otherwise )
        return Add(*args).as_leading_term(x)

    @classmethod
    def __eval_cond(cls, cond):
        """ Returns cond if it's a boolean, otherwise it is undecidable and returns None. """
        if isinstance(cond,bool):
            return cond
        return None

    @classmethod
    def _collapse_piecewise_args(cls, ecs, oth):
        """
        Collapse nested Piecewise expressions

        If an (expr,cond) pair in ecs is a Piecewise and the conditions of expr
        match the conditions in ecs, collapse expr to the expression matching
        cond.

        This allows the collapsing of Piecewise((Piecewise((x, x<0)), x<0)) to
        Piecewise((x, x<0)).

        Eventually, this code should be able to collapse Piecewise's having
        different intervals, but this will probably require using the new
        assumptions.
        """
        from sympy import Or
        # True when args are changed and we should rerun to check for new Piecewise exprs
        new_ecs = []
        or1 = Or( *[c for (_, c) in ecs] )
        for expr, cond in ecs:
            if isinstance(expr, Piecewise):
                or2 = Or( *[c for (_, c) in expr.exprcondpairs] )
                for e, c in expr.exprcondpairs:
                    # Don't collapse if cond is "True" as this leads to
                    # incorrect simplifications with nested Piecewises.
                    if c == cond and (or1 == or2 or cond is not True):
                        expr = e
            new_ecs.append( Tuple(expr, cond) )
        if isinstance(oth, Piecewise):
            or2 = Or( *[c for (_, c) in oth.exprcondpairs] )
            if or1 == or2:
                oth = oth.otherwise
        return new_ecs, oth


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
