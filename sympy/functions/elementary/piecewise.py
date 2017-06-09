from __future__ import print_function, division

from sympy.core import Basic, S, Function, diff, Tuple
from sympy.core.relational import Equality, Relational, _canonical
from sympy.functions.elementary.miscellaneous import Max, Min
from sympy.logic.boolalg import (And, Boolean, distribute_and_over_or, Not, Or,
    true, false)
from sympy.utilities.iterables import cartes
from sympy.core.compatibility import default_sort_key, range


class ExprCondPair(Tuple):
    """Represents an expression, condition pair."""

    def __new__(cls, expr, cond):
        if cond == True:
            return Tuple.__new__(cls, expr, true)
        elif cond == False:
            return Tuple.__new__(cls, expr, false)
        return Tuple.__new__(cls, expr, cond)

    @property
    def expr(self):
        """
        Returns the expression of this pair.
        """
        return self.args[0]

    @property
    def cond(self):
        """
        Returns the condition of this pair.
        """
        return self.args[1]

    @property
    def free_symbols(self):
        """
        Return the free symbols of this pair.
        """
        # Overload Basic.free_symbols because self.args[1] may contain non-Basic
        result = self.expr.free_symbols
        if hasattr(self.cond, 'free_symbols'):
            result |= self.cond.free_symbols
        return result

    @property
    def is_commutative(self):
        return self.expr.is_commutative

    def __iter__(self):
        yield self.expr
        yield self.cond


class Piecewise(Function):
    """
    Represents a piecewise function.

    Usage:

      Piecewise( (expr,cond), (expr,cond), ... )
        - Each argument is a 2-tuple defining an expression and condition
        - The conds are evaluated in turn returning the first that is True.
          If any of the evaluated conds are not determined explicitly False,
          e.g. x < 1, the function is returned in symbolic form.
        - If the function is evaluated at a place where all conditions are False,
          a ValueError exception will be raised.
        - Pairs where the cond is explicitly False, will be removed.

    Examples
    ========

      >>> from sympy import Piecewise, log
      >>> from sympy.abc import x
      >>> f = x**2
      >>> g = log(x)
      >>> p = Piecewise( (0, x<-1), (f, x<=1), (g, True))
      >>> p.subs(x,1)
      1
      >>> p.subs(x,5)
      log(5)

    See Also
    ========

    piecewise_fold
    """

    nargs = None
    is_Piecewise = True

    def __new__(cls, *args, **options):
        # (Try to) sympify args first
        newargs = []
        for ec in args:
            # ec could be a ExprCondPair or a tuple
            pair = ExprCondPair(*getattr(ec, 'args', ec))
            cond = pair.cond
            if cond == false:
                continue
            if not isinstance(cond, (bool, Relational, Boolean)):
                raise TypeError(
                    "Cond %s is of type %s, but must be a Relational,"
                    " Boolean, or a built-in bool." % (cond, type(cond)))
            newargs.append(pair)
            if cond == True:
                break

        if options.pop('evaluate', True):
            r = cls.eval(*newargs)
        else:
            r = None

        if r is None:
            return Basic.__new__(cls, *newargs, **options)
        else:
            return r

    @classmethod
    def eval(cls, *_args):
        """Either return a modified version of the args or, if no
        modifications were made, return None.

        Modifications that are made here:
        1) relationals are made canonical
        2) any False conditions are dropped
        3) any repeat of a previous condition is ignored
        3) any args past one with a true condition are dropped

        If there are no args left, an empty Piecewise will be returned.
        If there is a single arg with a True condition, its
        corresponding expression will be returned.
        """

        if not _args:
            return

        if len(_args) == 1 and _args[0][-1] == True:
            return _args[0][0]

        newargs = []  # the unevaluated conditions
        current_cond = set()  # the conditions up to a given e, c pair
        # make conditions canonical
        args = []
        for e, c in _args:
            if not c.is_Atom and not isinstance(c, Relational):
                free = c.free_symbols
                if len(free) == 1:
                    funcs = [i for i in c.atoms(Function)
                        if not isinstance(i, Boolean)]
                    if len(funcs) == 1 and len(
                            c.xreplace({list(funcs)[0]: Dummy()}
                            ).free_symbols) == 1:
                        # we can treat function like a symbol
                        free = funcs
                    _c = c
                    x = free.pop()
                    try:
                        c = c.as_set().as_relational(x)
                    except NotImplementedError:
                        pass
                    else:
                        reps = {}
                        for i in c.atoms(Relational):
                            ic = i.canonical
                            if ic.rhs in (S.Infinity, S.NegativeInfinity):
                                if not _c.has(ic.rhs):
                                    # don't accept introduction of
                                    # new Relationals with +/-oo
                                    reps[i] = S.true
                                elif ('=' not in ic.rel_op and
                                        c.xreplace({x: i.rhs}) !=
                                        _c.xreplace({x: i.rhs})):
                                    reps[i] = Relational(
                                        i.lhs, i.rhs, i.rel_op + '=')
                        c = c.xreplace(reps)
            args.append((e, _canonical(c)))

        for expr, cond in args:
            # Check here if expr is a Piecewise and collapse if one of
            # the conds in expr matches cond. This allows the collapsing
            # of Piecewise((Piecewise((x,x<0)),x<0)) to Piecewise((x,x<0)).
            # This is important when using piecewise_fold to simplify
            # multiple Piecewise instances having the same conds.
            # Eventually, this code should be able to collapse Piecewise's
            # having different intervals, but this will probably require
            # using the new assumptions.
            if isinstance(expr, Piecewise):
                unmatching = []
                for i, (e, c) in enumerate(expr.args):
                    if c in current_cond:
                        # this would already have triggered
                        continue
                    if c == cond:
                        if c != True:
                            # nothing past this condition will ever
                            # trigger and only those args before this
                            # that didn't match a previous condition
                            # could possibly trigger
                            if unmatching:
                                expr = Piecewise(*(
                                    unmatching + [(e, c)]))
                            else:
                                expr = e
                        break
                    else:
                        unmatching.append((e, c))

            # check for condition repeats
            got = False
            # -- if an And contains a condition that was
            #    already encountered, then the And will be
            #    False: if the previous condition was False
            #    then the And will be False and if the previous
            #    condition is True then then we wouldn't get to
            #    this point. In either case, we can skip this condition.
            for i in ([cond] +
                    (list(cond.args) if isinstance(cond, And) else
                    [])):
                if i in current_cond:
                    got = True
                    break
            if got:
                continue

            # -- if not(c) is already in current_cond then c is
            #    a redundant condition in an And. This does not
            #    apply to Or, however: (e1, c), (e2, Or(~c, d))
            #    is not (e1, c), (e2, d) because if c and d are
            #    both False this would give no results when the
            #    true answer should be (e2, True)
            if isinstance(cond, And):
                nonredundant = []
                for c in cond.args:
                    if (isinstance(c, Relational) and
                            (~c).canonical in current_cond):
                        continue
                    nonredundant.append(c)
                cond = cond.func(*nonredundant)
            elif isinstance(cond, Relational):
                if (~cond).canonical in current_cond:
                    cond = S.true

            current_cond.add(cond)

            # collect successive e,c pairs when exprs or cond match
            if newargs:
                if newargs[-1].expr == expr:
                    orcond = Or(cond, newargs[-1].cond)
                    if isinstance(orcond, (And, Or)):
                        orcond = distribute_and_over_or(orcond)
                    newargs[-1] = ExprCondPair(expr, orcond)
                    continue
                elif newargs[-1].cond == cond:
                    orexpr = Or(expr, newargs[-1].expr)
                    if isinstance(orexpr, (And, Or)):
                        orexpr = distribute_and_over_or(orexpr)
                    newargs[-1] == ExprCondPair(orexpr, cond)
                    continue

            newargs.append(ExprCondPair(expr, cond))

        # some conditions may have been redundant
        missing = len(newargs) != len(_args)
        # some conditions may have changed
        same = all(a == b for a, b in zip(newargs, _args))
        # if either change happened we return the expr with the
        # updated args
        if not newargs:
            raise ValueError(filldedent('''
There are no conditions (or none that are not trivially false) to
define an expression.'''))
        if missing or not same:
            return cls(*newargs)

    def doit(self, **hints):
        """
        Evaluate this piecewise function.
        """
        newargs = []
        for e, c in self.args:
            if hints.get('deep', True):
                if isinstance(e, Basic):
                    e = e.doit(**hints)
                if isinstance(c, Basic):
                    c = c.doit(**hints)
            newargs.append((e, c))
        return self.func(*newargs)

    def _eval_as_leading_term(self, x):
        for e, c in self.args:
            if c == True or c.subs(x, 0) == True:
                return e.as_leading_term(x)

    def _eval_adjoint(self):
        return self.func(*[(e.adjoint(), c) for e, c in self.args])

    def _eval_conjugate(self):
        return self.func(*[(e.conjugate(), c) for e, c in self.args])

    def _eval_derivative(self, x):
        return self.func(*[(diff(e, x), c) for e, c in self.args])

    def _eval_evalf(self, prec):
        return self.func(*[(e.evalf(prec), c) for e, c in self.args])

    def _eval_integral(self, x):
        from sympy.integrals import integrate
        return self.func(*[(integrate(e, x), c) for e, c in self.args])

    def _eval_interval(self, sym, a, b):
        """Evaluates the function along the sym in a given interval ab"""
        # FIXME: Currently complex intervals are not supported.  A possible
        # replacement algorithm, discussed in issue 5227, can be found in the
        # following papers;
        #     http://portal.acm.org/citation.cfm?id=281649
        #     http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.70.4127&rep=rep1&type=pdf

        if a is None or b is None:
            # In this case, it is just simple substitution
            return piecewise_fold(
                super(Piecewise, self)._eval_interval(sym, a, b))

        mul = 1
        if (a == b) == True:
            return S.Zero
        elif (a > b) == True:
            a, b, mul = b, a, -1
        elif (a <= b) != True:
            newargs = []
            for e, c in self.args:
                intervals = self._sort_expr_cond(
                    sym, S.NegativeInfinity, S.Infinity, c)
                values = []
                for lower, upper, expr in intervals:
                    if (a < lower) == True:
                        mid = lower
                        rep = b
                        val = e._eval_interval(sym, mid, b)
                        val += self._eval_interval(sym, a, mid)
                    elif (a > upper) == True:
                        mid = upper
                        rep = b
                        val = e._eval_interval(sym, mid, b)
                        val += self._eval_interval(sym, a, mid)
                    elif (a >= lower) == True and (a <= upper) == True:
                        rep = b
                        val = e._eval_interval(sym, a, b)
                    elif (b < lower) == True:
                        mid = lower
                        rep = a
                        val = e._eval_interval(sym, a, mid)
                        val += self._eval_interval(sym, mid, b)
                    elif (b > upper) == True:
                        mid = upper
                        rep = a
                        val = e._eval_interval(sym, a, mid)
                        val += self._eval_interval(sym, mid, b)
                    elif ((b >= lower) == True) and ((b <= upper) == True):
                        rep = a
                        val = e._eval_interval(sym, a, b)
                    else:
                        raise NotImplementedError(
                            """The evaluation of a Piecewise interval when both the lower
                            and the upper limit are symbolic is not yet implemented.""")
                    values.append(val)
                if len(set(values)) == 1:
                    try:
                        c = c.subs(sym, rep)
                    except AttributeError:
                        pass
                    e = values[0]
                    newargs.append((e, c))
                else:
                    for i in range(len(values)):
                        newargs.append((values[i], (c == True and i == len(values) - 1) or
                            And(rep >= intervals[i][0], rep <= intervals[i][1])))
            return self.func(*newargs)

        # Determine what intervals the expr,cond pairs affect.
        int_expr = self._sort_expr_cond(sym, a, b)

        # Finally run through the intervals and sum the evaluation.
        ret_fun = 0
        for int_a, int_b, expr in int_expr:
            if isinstance(expr, Piecewise):
                # If we still have a Piecewise by now, _sort_expr_cond would
                # already have determined that its conditions are independent
                # of the integration variable, thus we just use substitution.
                ret_fun += piecewise_fold(
                    super(Piecewise, expr)._eval_interval(sym, Max(a, int_a), Min(b, int_b)))
            else:
                ret_fun += expr._eval_interval(sym, Max(a, int_a), Min(b, int_b))
        return mul * ret_fun

    def _sort_expr_cond(self, sym, a, b, targetcond=None):
        """Determine what intervals the expr, cond pairs affect.

        1) If cond is True, then log it as default
        1.1) Currently if cond can't be evaluated, throw NotImplementedError.
        2) For each inequality, if previous cond defines part of the interval
           update the new conds interval.
           -  eg x < 1, x < 3 -> [oo,1],[1,3] instead of [oo,1],[oo,3]
        3) Sort the intervals to make it easier to find correct exprs

        Under normal use, we return the expr,cond pairs in increasing order
        along the real axis corresponding to the symbol sym.  If targetcond
        is given, we return a list of (lowerbound, upperbound) pairs for
        this condition."""
        from sympy.solvers.inequalities import _solve_inequality
        default = None
        int_expr = []
        expr_cond = []
        or_cond = False
        or_intervals = []
        independent_expr_cond = []
        for expr, cond in self.args:
            if isinstance(cond, Or):
                for cond2 in sorted(cond.args, key=default_sort_key):
                    expr_cond.append((expr, cond2))
            else:
                expr_cond.append((expr, cond))
            if cond == True:
                break
        for expr, cond in expr_cond:
            if cond == True:
                independent_expr_cond.append((expr, cond))
                default = self.func(*independent_expr_cond)
                break
            orig_cond = cond
            if sym not in cond.free_symbols:
                independent_expr_cond.append((expr, cond))
                continue
            elif isinstance(cond, Equality):
                continue
            elif isinstance(cond, And):
                lower = S.NegativeInfinity
                upper = S.Infinity
                for cond2 in cond.args:
                    if sym not in [cond2.lts, cond2.gts]:
                        cond2 = _solve_inequality(cond2, sym)
                    if cond2.lts == sym:
                        upper = Min(cond2.gts, upper)
                    elif cond2.gts == sym:
                        lower = Max(cond2.lts, lower)
                    else:
                        raise NotImplementedError(
                            "Unable to handle interval evaluation of expression.")
            else:
                if sym not in [cond.lts, cond.gts]:
                    cond = _solve_inequality(cond, sym)
                lower, upper = cond.lts, cond.gts  # part 1: initialize with givens
                if cond.lts == sym:                # part 1a: expand the side ...
                    lower = S.NegativeInfinity   # e.g. x <= 0 ---> -oo <= 0
                elif cond.gts == sym:            # part 1a: ... that can be expanded
                    upper = S.Infinity           # e.g. x >= 0 --->  oo >= 0
                else:
                    raise NotImplementedError(
                        "Unable to handle interval evaluation of expression.")

            # part 1b: Reduce (-)infinity to what was passed in.
            lower, upper = Max(a, lower), Min(b, upper)

            for n in range(len(int_expr)):
                # Part 2: remove any interval overlap.  For any conflicts, the
                # iterval already there wins, and the incoming interval updates
                # its bounds accordingly.
                if self.__eval_cond(lower < int_expr[n][1]) and \
                        self.__eval_cond(lower >= int_expr[n][0]):
                    lower = int_expr[n][1]
                elif len(int_expr[n][1].free_symbols) and \
                        self.__eval_cond(lower >= int_expr[n][0]):
                    if self.__eval_cond(lower == int_expr[n][0]):
                        lower = int_expr[n][1]
                    else:
                        int_expr[n][1] = Min(lower, int_expr[n][1])
                elif len(int_expr[n][0].free_symbols) and \
                        self.__eval_cond(upper == int_expr[n][1]):
                    upper = Min(upper, int_expr[n][0])
                elif len(int_expr[n][1].free_symbols) and \
                        (lower >= int_expr[n][0]) != True and \
                        (int_expr[n][1] == Min(lower, upper)) != True:
                    upper = Min(upper, int_expr[n][0])
                elif self.__eval_cond(upper > int_expr[n][0]) and \
                        self.__eval_cond(upper <= int_expr[n][1]):
                    upper = int_expr[n][0]
                elif len(int_expr[n][0].free_symbols) and \
                        self.__eval_cond(upper < int_expr[n][1]):
                    int_expr[n][0] = Max(upper, int_expr[n][0])

            if self.__eval_cond(lower >= upper) != True:  # Is it still an interval?
                int_expr.append([lower, upper, expr])
            if orig_cond == targetcond:
                return [(lower, upper, None)]
            elif isinstance(targetcond, Or) and cond in targetcond.args:
                or_cond = Or(or_cond, cond)
                or_intervals.append((lower, upper, None))
                if or_cond == targetcond:
                    or_intervals.sort(key=lambda x: x[0])
                    return or_intervals

        int_expr.sort(key=lambda x: x[1].sort_key(
        ) if x[1].is_number else S.NegativeInfinity.sort_key())
        int_expr.sort(key=lambda x: x[0].sort_key(
        ) if x[0].is_number else S.Infinity.sort_key())

        for n in range(len(int_expr)):
            if len(int_expr[n][0].free_symbols) or len(int_expr[n][1].free_symbols):
                if isinstance(int_expr[n][1], Min) or int_expr[n][1] == b:
                    newval = Min(*int_expr[n][:-1])
                    if n > 0 and int_expr[n][0] == int_expr[n - 1][1]:
                        int_expr[n - 1][1] = newval
                    int_expr[n][0] = newval
                else:
                    newval = Max(*int_expr[n][:-1])
                    if n < len(int_expr) - 1 and int_expr[n][1] == int_expr[n + 1][0]:
                        int_expr[n + 1][0] = newval
                    int_expr[n][1] = newval

        # Add holes to list of intervals if there is a default value,
        # otherwise raise a ValueError.
        holes = []
        curr_low = a
        for int_a, int_b, expr in int_expr:
            if (curr_low < int_a) == True:
                holes.append([curr_low, Min(b, int_a), default])
            elif (curr_low >= int_a) != True:
                holes.append([curr_low, Min(b, int_a), default])
            curr_low = Min(b, int_b)
        if (curr_low < b) == True:
            holes.append([Min(b, curr_low), b, default])
        elif (curr_low >= b) != True:
            holes.append([Min(b, curr_low), b, default])

        if holes and default is not None:
            int_expr.extend(holes)
            if targetcond == True:
                return [(h[0], h[1], None) for h in holes]
        elif holes and default is None:
            raise ValueError("Called interval evaluation over piecewise "
                             "function on undefined intervals %s" %
                             ", ".join([str((h[0], h[1])) for h in holes]))

        return int_expr

    def _eval_nseries(self, x, n, logx):
        args = [(ec.expr._eval_nseries(x, n, logx), ec.cond) for ec in self.args]
        return self.func(*args)

    def _eval_power(self, s):
        return self.func(*[(e**s, c) for e, c in self.args])

    def _eval_subs(self, old, new):
        """
        Piecewise conditions may contain bool which are not of Basic type.
        """
        args = list(self.args)
        for i, (e, c) in enumerate(args):
            if isinstance(c, bool):
                pass
            elif isinstance(c, Basic):
                c = c._subs(old, new)
            if c != False:
                e = e._subs(old, new)
            args[i] = e, c
            if c == True:
                return self.func(*args)

        return self.func(*args)

    def _eval_transpose(self):
        return self.func(*[(e.transpose(), c) for e, c in self.args])

    def _eval_template_is_attr(self, is_attr, when_multiple=None):
        b = None
        for expr, _ in self.args:
            a = getattr(expr, is_attr)
            if a is None:
                return None
            if b is None:
                b = a
            elif b is not a:
                return when_multiple
        return b

    _eval_is_finite = lambda self: self._eval_template_is_attr(
        'is_finite', when_multiple=False)
    _eval_is_complex = lambda self: self._eval_template_is_attr('is_complex')
    _eval_is_even = lambda self: self._eval_template_is_attr('is_even')
    _eval_is_imaginary = lambda self: self._eval_template_is_attr(
        'is_imaginary')
    _eval_is_integer = lambda self: self._eval_template_is_attr('is_integer')
    _eval_is_irrational = lambda self: self._eval_template_is_attr(
        'is_irrational')
    _eval_is_negative = lambda self: self._eval_template_is_attr('is_negative')
    _eval_is_nonnegative = lambda self: self._eval_template_is_attr(
        'is_nonnegative')
    _eval_is_nonpositive = lambda self: self._eval_template_is_attr(
        'is_nonpositive')
    _eval_is_nonzero = lambda self: self._eval_template_is_attr(
        'is_nonzero', when_multiple=True)
    _eval_is_odd = lambda self: self._eval_template_is_attr('is_odd')
    _eval_is_polar = lambda self: self._eval_template_is_attr('is_polar')
    _eval_is_positive = lambda self: self._eval_template_is_attr('is_positive')
    _eval_is_real = lambda self: self._eval_template_is_attr('is_real')
    _eval_is_zero = lambda self: self._eval_template_is_attr(
        'is_zero', when_multiple=False)

    @classmethod
    def __eval_cond(cls, cond):
        """Return the truth value of the condition."""
        from sympy.solvers.solvers import checksol
        if cond == True:
            return True
        if isinstance(cond, Equality):
            diff = cond.lhs - cond.rhs
            if diff.is_commutative:
                return diff.is_zero
        return None

    def as_expr_set_pairs(self):
        exp_sets = []
        U = S.Reals
        for expr, cond in self.args:
            cond_int = U.intersect(cond.as_set())
            U = U - cond_int
            exp_sets.append((expr, cond_int))
        return exp_sets


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
    Piecewise((x**2, x < 1), (x, True))

    See Also
    ========

    Piecewise
    """
    if not isinstance(expr, Basic) or not expr.has(Piecewise):
        return expr
    if isinstance(expr, (ExprCondPair, Piecewise)):
        new_args = []
        for e, c in expr.args:
            if not isinstance(e, Piecewise):
                e = piecewise_fold(e)
            if isinstance(e, Piecewise):
                new_args.extend([(piecewise_fold(ei), And(ci, c)) for ei, ci in e.args])
            else:
                new_args.append((e, c))
    else:
        new_args = []
        folded = list(map(piecewise_fold, expr.args))
        for ec in cartes(*[
                (i.args if isinstance(i, Piecewise) else
                 [(i, S.true)]) for i in folded]):
            e, c = zip(*ec)
            new_args.append((expr.func(*e), And(*c)))
    rv = Piecewise(*new_args)
    return rv
