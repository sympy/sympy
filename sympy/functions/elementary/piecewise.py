from __future__ import print_function, division

from sympy.core import Basic, S, Function, diff, Tuple, Dummy
from sympy.core.sympify import _sympify
from sympy.core.relational import Equality, Relational, _canonical
from sympy.functions.elementary.miscellaneous import Max, Min
from sympy.logic.boolalg import (And, Boolean, distribute_and_over_or, Not, Or,
    true, false)
from sympy.utilities.iterables import cartes
from sympy.core.compatibility import default_sort_key, range
from sympy.utilities.iterables import uniq, is_sequence, ordered, product
from sympy.utilities.misc import filldedent, Undecidable


Undefined = S.NaN  # Piecewise()

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
        return getattr(self.expr, 'is_commutative', False)

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

    def piecewise_integrate(self, x, **kwargs):
        """Return the Piecewise with each expression being
        replaced with its antiderivative. To obtain a continuous
        antiderivative, use the `integrate` function or method.

        Examples
        ========

        >>> from sympy import Piecewise
        >>> from sympy.abc import x
        >>> p = Piecewise((0, x < 0), (1, x < 1), (2, True))
        >>> p.piecewise_integrate(x)
        Piecewise((0, x < 0), (x, x < 1), (2*x, True))

        Note that this does not give a continuous function, e.g.
        at x = 1 the 3rd condition applies and the antiderivative
        there is 2*x so the value of the antiderivative is 2:

        >>> anti = _
        >>> anti.subs(x, 1)
        2

        The continuous derivative accounts for the integral *up to*
        the point of interest, however:

        >>> p.integrate(x)
        Piecewise((0, x < 0), (x, x < 1), (2*x - 1, True))
        >>> _.subs(x, 1)
        1

        See Also
        ========
        Piecewise._eval_integral
        """
        from sympy.integrals import integrate
        return self.func(*[(integrate(e, x, **kwargs), c) for e, c in self.args])

    def _handle_irel(self, x, handler):
        """Return either None (if the conditions of self depend only on x) else
        a Piecewise expression whose expressions (handled by the handler that
        was passed) are paired with the governing x-independent relationals,
        e.g. Piecewise((A, a(x) & b(y)), (B, c(x) | c(y)) ->
        Piecewise(
            (handler(Piecewise((A, a(x) & True), (B, c(x) | True)), b(y) & c(y)),
            (handler(Piecewise((A, a(x) & True), (B, c(x) | False)), b(y)),
            (handler(Piecewise((A, a(x) & False), (B, c(x) | True)), c(y)),
            (handler(Piecewise((A, a(x) & False), (B, c(x) | False)), True))
        """
        # identify governing relationals
        rel = self.atoms(Relational)
        irel = list(ordered([r for r in rel if x not in r.free_symbols
            and r not in (S.true, S.false)]))
        if irel:
            args = {}
            exprinorder = []
            for truth in product((1, 0), repeat=len(irel)):
                reps = dict(zip(irel, truth))
                # only store the true conditions since the false are implied
                # when they appear lower in the Piecewise args
                if 1 not in truth:
                    cond = None  # flag this one so it doesn't get combined
                else:
                    andargs = Tuple(*[i for i in reps if reps[i]])
                    free = list(andargs.free_symbols)
                    if len(free) == 1:
                        from sympy.solvers.inequalities import reduce_inequalities
                        t = reduce_inequalities(andargs, free[0])
                    else:
                        t = And(*andargs)
                    if t is S.false:
                        continue  # an impossible combination
                    cond = t
                expr = handler(self.xreplace(reps))
                if isinstance(expr, self.func) and len(expr.args) == 1:
                    expr, econd = expr.args[0]
                    cond = And(econd, True if cond is None else cond)
                # the ec pairs are being collected since all possibilities
                # are being enumerated, but don't put the last one in since
                # its expr might match a previous expression and it
                # must appear last in the args
                if cond is not None:
                    args.setdefault(expr, []).append(cond)
                    # but since we only store the true conditions we must maintain
                    # the order so that the expression with the most true values
                    # comes first
                    exprinorder.append(expr)
            # convert collected conditions as args of Or
            for k in args:
                args[k] = Or(*args[k])
            # take them in the order obtained
            args = [(e, args[e]) for e in uniq(exprinorder)]
            # add in the last arg
            args.append((expr, True))
            # if any condition reduced to True, it needs to go last
            # and there should only be one of them or else the exprs
            # should agree
            trues = [i for i in range(len(args)) if args[i][1] is S.true]
            if not trues:
                # make the last one True since all cases were enumerated
                e, c = args[-1]
                args[-1] = (e, S.true)
            else:
                assert len(set([e for e, c in [args[i] for i in trues]])) == 1
                args.append(args.pop(trues.pop()))
                while trues:
                    args.pop(trues.pop())
            return Piecewise(*args)

    def _eval_integral(self, x, _first=True, **kwargs):
        """Return the indefinite integral of the
        Piecewise such that subsequent substitution of x with a
        value will give the value of the integral (not including
        the constant of integration) up to that point. To only
        integrate the individual parts of Piecewise, use the
        `piecewise_integrate` method.

        Examples
        ========

        >>> from sympy import Piecewise
        >>> from sympy.abc import x
        >>> p = Piecewise((0, x < 0), (1, x < 1), (2, True))
        >>> p.integrate(x)
        Piecewise((0, x < 0), (x, x < 1), (2*x - 1, True))
        >>> p.piecewise_integrate(x)
        Piecewise((0, x < 0), (x, x < 1), (2*x, True))

        See Also
        ========
        Piecewise.piecewise_integrate
        """
        from sympy.integrals.integrals import integrate

        if _first:
            def handler(ipw):
                if isinstance(ipw, self.func):
                    return ipw._eval_integral(x, _first=False, **kwargs)
                else:
                    return ipw.integrate(x, **kwargs)
            irv = self._handle_irel(x, handler)
            if irv is not None:
                return irv

        # solve a Piecewise with lo <= hi and no x-independent relationals
        # ------------------------------------------------------------------
        try:
            abei = self._abei(x)
        except NotImplementedError:
            from sympy import Integral
            return Integral(self, (x, lo, hi))  # unevaluated

        pieces = [(a, b) for a, b, _, _ in abei]
        oo = S.Infinity
        done = [(-oo, oo, -1)]
        for k, p in enumerate(pieces):
            if p[:2] == (-oo, oo):
                # all undone intervals will get this key
                for j, (a, b, i) in enumerate(done):
                    if i == -1:
                        done[j] = a, b, k
                break  # nothing else to consider
            N = len(done) - 1
            for j, (a, b, i) in enumerate(reversed(done)):
                if i == -1:
                    j = N - j
                    done[j: j + 1] = _clip(p, (a, b), k)

        # check for holes
        for a, b, i in done:
            if i == -1 and a != b:
                # when we reference arg -1 we will get Undefined
                abei.append((-oo, oo, Undefined, -1))
                break

        # return the sum of the intervals
        args = []
        sum = None
        for a, b, i in done:
            anti = integrate(abei[i][-2], x, **kwargs)
            if sum is None:
                sum = anti
            else:
                sum = sum.subs(x, a)
                sum += anti._eval_interval(x, a, x)
            # see if we know whether b is contained in original
            # condition
            if b is S.Infinity:
                cond = True
            elif self.args[abei[i][-1]].cond.subs(x, b) == False:
                cond = (x < b)
            else:
                cond = (x <= b)
            args.append((_mmsimp(sum), cond))
        return Piecewise(*args)

    def _eval_interval(self, sym, a, b, _first=True):
        """Evaluates the function along the sym in a given interval [a, b]"""
        # FIXME: Currently complex intervals are not supported.  A possible
        # replacement algorithm, discussed in issue 5227, can be found in the
        # following papers;
        #     http://portal.acm.org/citation.cfm?id=281649
        #     http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.70.4127&rep=rep1&type=pdf
        from sympy.core.symbol import Dummy

        if a is None or b is None:
            # In this case, it is just simple substitution
            return super(Piecewise, self)._eval_interval(sym, a, b)
        else:
            x, lo, hi = map(_sympify, (sym, a, b))

        if _first:  # get only x-dependent relationals
            def handler(ipw):
                if isinstance(ipw, self.func):
                    return ipw._eval_interval(x, lo, hi, _first=None)
                else:
                    return ipw._eval_interval(x, lo, hi)
            irv = self._handle_irel(x, handler)
            if irv is not None:
                return irv

            if (lo < hi) is S.true or hi is S.Infinity or lo is S.NegativeInfinity:
                pass
            else:
                neg = -self._eval_interval(x, hi, lo, _first=False)
                if (lo < hi) is S.false or lo is S.Infinity or hi is S.NegativeInfinity:
                    rv = neg
                else:
                    rv = Piecewise(
                        (self._eval_interval(x, lo, hi, _first=False),
                            lo < hi),
                        (neg,
                            True))
                if rv == Undefined:
                    raise ValueError("Can't integrate across undefined region.")
                return rv

        # solve a Piecewise with lo <= hi and no x-independent relationals
        # ------------------------------------------------------------------
        try:
            abei = self._abei(x)
        except NotImplementedError:
            from sympy import Integral
            return Integral(self, (x, lo, hi))  # unevaluated

        pieces = [(a, b) for a, b, _, _ in abei]
        done = [(lo, hi, -1)]
        oo = S.Infinity
        for k, p in enumerate(pieces):
            if p[:2] == (-oo, oo):
                # all undone intervals will get this key
                for j, (a, b, i) in enumerate(done):
                    if i == -1:
                        done[j] = a, b, k
                break  # nothing else to consider
            N = len(done) - 1
            for j, (a, b, i) in enumerate(reversed(done)):
                if i == -1:
                    j = N - j
                    done[j: j + 1] = _clip(p, (a, b), k)
        # check for holes
        for a, b, i in done:
            if i == -1:
                if a != b:
                    # this requires a default (i == -1)
                    # over a potentially non-zero interval
                    # since a != b so we can't do this
                    return S.NaN
                else:
                    # XXX are there expressions that would give
                    # a value over a zero-width interval?
                    # Integrating DiracDelta or Heaviside
                    # does but we are not integrating right
                    # now. If there are things to check for
                    # check expression given by abei[i][-2]
                    pass
        # return the sum of the intervals
        sum = S.Zero
        for a, b, i in done:
            sum += abei[i][-2]._eval_interval(x, a, b)
        return _mmsimp(sum)

    def _abei(self, sym):
        """Return (a, b, e, i) where a and b are the lower and
        upper bounds in which the expression e of argument i
        in self is defined.

        If there are any relationals not involving sym, or
        any relational cannot be solved for sym, we raise
        NotImplementedError. The evaluated conditions will
        be returned as ranges. Discontinuous ranges will be
        returned separately with identical expressions and
        the first condition that evaluates to True will be
        returned as the last tuple with a, b = -oo, oo.
        """
        from sympy.solvers.inequalities import _solve_inequality
        from sympy.logic.boolalg import to_cnf, distribute_or_over_and

        assert isinstance(self, Piecewise)

        def _solve_relational(r):
            rv = _solve_inequality(r, sym)
            if isinstance(rv, Relational) and \
                    sym in rv.free_symbols:
                if rv.args[0] != sym:
                    raise NotImplementedError(filldedent('''
Unable to solve relational %s for %s.''' % (r, sym)))
                if rv.rel_op == '!=':
                    rv = Or(sym < rv.rhs, sym > rv.rhs)
            if rv == (S.NegativeInfinity < sym) & (sym < S.Infinity):
                rv = S.true
            return rv

        def nonsymfail(cond):
                    raise NotImplementedError(filldedent('''
A condition not involving %s appeared: %s''' % (sym, cond)))

        # make self canonical wrt Relationals
        reps = dict(
            [(r,_solve_relational(r)) for r in self.atoms(Relational)])
        # process args individually so if any evaluate, their position
        # in the original Piecewise will be known
        args = [i.xreplace(reps) for i in self.args]

        # precondition args
        expr_cond = []
        default = idefault = None
        for i, (expr, cond) in enumerate(args):
            if cond == False:
                continue
            elif cond == True:
                default = expr
                idefault = i
                break
            elif sym not in cond.free_symbols:
                nonsymfail(cond)

            cond = to_cnf(cond)
            if isinstance(cond, And):
                cond = distribute_or_over_and(cond)

            if isinstance(cond, Or):
                expr_cond.extend(
                    [(i, expr, o) for o in cond.args])
            elif cond != False:
                expr_cond.append((i, expr, cond))

        # determine intervals represented by conditions
        int_expr = []
        for iarg, expr, cond in expr_cond:
            if isinstance(cond, Equality):
                if cond.lhs == sym:
                    v = cond.rhs
                    int_expr.append((v, v, expr.subs(sym, v), iarg))
                    continue
                else:
                    nonsymfail(cond)

            if isinstance(cond, And):
                lower = S.NegativeInfinity
                upper = S.Infinity
                nonsym = False
                rhs = None
                for cond2 in cond.args:
                    if isinstance(cond2, Equality):
                        # all relationals were solved and
                        # only sym-dependent ones made it to here
                        assert cond2.lhs == sym
                        assert sym not in cond2.rhs.free_symbols
                        if rhs is None:
                            rhs = cond2.rhs
                        else:
                            same = Eq(cond2.rhs, rhs)
                            if same == False:
                                cond = S.false
                                break
                            elif same != True:
                                raise Undecidable('%s did not evaluate' % same)
                    elif cond2.lts == sym:
                        upper = Min(cond2.gts, upper)
                    elif cond2.gts == sym:
                        lower = Max(cond2.lts, lower)
                    else:
                        # mark but don't fail until later
                        nonsym = cond2
                if rhs is not None and cond not in (S.true, S.false):
                    lower = upper = rhs
                    cond = cond.subs(sym, lower)
                # cond might have evaluated b/c of an Eq
                if cond is S.true:
                    default = expr
                    continue
                elif cond is S.false:
                    continue
                elif nonsym is not False:
                    nonsymfail(nonsym)
                else:
                    pass  # for coverage
            elif isinstance(cond, Relational):
                lower, upper = cond.lts, cond.gts  # part 1: initialize with givens
                if cond.lts == sym:                # part 1a: expand the side ...
                    lower = S.NegativeInfinity   # e.g. x <= 0 ---> -oo <= 0
                elif cond.gts == sym:            # part 1a: ... that can be expanded
                    upper = S.Infinity           # e.g. x >= 0 --->  oo >= 0
                else:
                    nonsymfail(cond)
            else:
                raise NotImplementedError(
                    'unrecognized condition: %s' % cond)

            lower, upper = lower, Max(lower, upper)
            if (lower > upper) is not S.true:
                int_expr.append((lower, upper, expr, iarg))

        if default is not None:
            int_expr.append(
                (S.NegativeInfinity, S.Infinity, default, idefault))

        return int_expr

    def _eval_nseries(self, x, n, logx):
        args = [(ec.expr._eval_nseries(x, n, logx), ec.cond) for ec in self.args]
        return self.func(*args)

    def _eval_power(self, s):
        return self.func(*[(e**s, c) for e, c in self.args])

    def _eval_subs(self, old, new):
        args = list(self.args)
        for i, (e, c) in enumerate(args):
            if isinstance(c, Basic):
                c = c._subs(old, new)
            if c != False:
                try:
                    e = e._subs(old, new)
                except AttributeError:
                    if is_sequence(e):
                        e = type(e)([ei.subs(old, new) for ei in e])
                    else:
                        pass
            args[i] = (e, c)
            if c == True:
                break

        return self.func(*args)

    def _eval_transpose(self):
        return self.func(*[(e.transpose(), c) for e, c in self.args])

    def _eval_template_is_attr(self, is_attr):
        b = None
        for expr, _ in self.args:
            a = getattr(expr, is_attr)
            if a is None:
                return
            if b is None:
                b = a
            elif b is not a:
                return
        return b

    _eval_is_finite = lambda self: self._eval_template_is_attr(
        'is_finite')
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
        'is_nonzero')
    _eval_is_odd = lambda self: self._eval_template_is_attr('is_odd')
    _eval_is_polar = lambda self: self._eval_template_is_attr('is_polar')
    _eval_is_positive = lambda self: self._eval_template_is_attr('is_positive')
    _eval_is_real = lambda self: self._eval_template_is_attr('is_real')
    _eval_is_zero = lambda self: self._eval_template_is_attr(
        'is_zero')

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

    new_args = []
    if isinstance(expr, (ExprCondPair, Piecewise)):
        for e, c in expr.args:
            if not isinstance(e, Piecewise):
                e = piecewise_fold(e)
            if isinstance(e, Piecewise):
                new_args.extend([(piecewise_fold(ei), And(ci, c))
                    for ei, ci in e.args])
            else:
                new_args.append((e, c))
    else:
        from sympy.utilities.iterables import cartes
        folded = list(map(piecewise_fold, expr.args))
        for ec in cartes(*[
                (i.args if isinstance(i, Piecewise) else
                 [(i, S.true)]) for i in folded]):
            e, c = zip(*ec)
            new_args.append((expr.func(*e), And(*c)))

    return Piecewise(*new_args)


def _pair(A, B, union=False, dict=None, abi=False):
    """Return the intervals on the real number line that are covered
    by the two intervals, A = [a, b] and B = [c, d]. If union is
    true then a simple list of intervals is returned; if dict is
    true (default) then the intervals are keyed to 0 if they are
    covered by A, 1 if they are covered by B and -1 if they are covered
    by neither, with higher precedence given to A. If abi is true,
    return an ordered set of tuples (a, b, i) where
    a is the lower bound, b is the upper bound and i indicates which
    interval covers that range: 0 for A, 1 for B and -1 for neither.

    Examples
    ========

    >>> from sympy import Tuple
    >>> from sympy.functions.elementary.piecewise import _pair
    >>> from sympy.abc import x

    For the following segments,

           1---2
                      4---5

    5 regions are identified:

        -oo -> 1 -> 2 -> 4 -> 5 -> oo

    >>> A = (1, 2); B = (4, 5)
    >>> _pair(A, B, abi=True)
    [(-oo, 1, -1), (1, 2, 0), (2, 4, -1), (4, 5, 1), (5, oo, -1)]
    >>> _pair(A, B, union=True)
    [(1, 2), (4, 5)]
    >>> _pair(A, B, dict=True)  # default
    {-1: [(-oo, 1), (2, 4), (5, oo)], 0: [(1, 2)], 1: [(4, 5)]}

    When one or more of the boundaries are symbolic, the results
    become more complex and the boundaries are given in terms of
    Min and Max functions:

    >>> _pair((x, 2), (4, 5))
    {-1: [(-oo, Min(4, x)), (Min(5, x), x), (2, 4), (5, oo)],
    0: [(x, 2)],
    1: [(Min(4, x), Min(5, x)), (4, 5)]}

    When an actual value is given for x, there may be redundant
    intervals and intervals with no width:

    >>> last = _
    >>> def replace(x, y):
    ...     rv = {}
    ...     for k, v in last.items():
    ...         rv[k] = [Tuple(*i).subs(x, y) for i in v]
    ...     return rv
    ...
    >>> for k, v in sorted(replace(x, 1).items()):
    ...     print(k, v)
    (-1, [(-oo, 1), (1, 1), (2, 4), (5, oo)])
    (0, [(1, 2)])
    (1, [(1, 1), (4, 5)])

    There may also be intervals wherein the right boundary is greater
    than the left and redundant intervals:

    >>> for k, v in sorted(replace(x, 6).items()):
    ...     print(k, v)
    (-1, [(-oo, 4), (5, 6), (2, 4), (5, oo)])
    (0, [(6, 2)])
    (1, [(4, 5), (4, 5)])
    """
    from sympy.logic.boolalg import Xor
    if not any(i for i in (abi, dict, union)):
        dict = True
    if [bool(i) for i in (abi, dict, union)].count(True) != 1:
        raise ValueError('only 1 of abi, dict and union may be True')

    a, b = A
    c, d = B
    oo = S.Infinity
    # the real line is broken into pieces which are compactly
    # represented as -oo, i1, hi1_lo2, i2, hi2_lo3, i3, ..., oo
    # where i1, i2, ...., refer to the interval A (0), B(1) or
    # neither (-1) as covering that intrval
    if a is -oo and b is oo:
        br = a, 0, b
    elif c is -oo and d is oo:
        if a is -oo:
            br = a, 0, b, 1, oo
        elif b is oo:
            br = -oo, 1, a, 0, b
        else:
            br = -oo, 1, a, 0, b, 1, oo
    elif b is oo and d is oo:  # (a, oo) (c, oo)
        br = -oo, -1, Min(a, c), 1, a, 0, oo
    elif a is -oo and c is -oo:  # (-oo, b) (-oo, d)
        br = -oo, 0, b, 1, Max(b, d), -1, oo
    elif a is -oo and d is oo:  # (-oo, b) (c, oo)
        br = -oo, 0, b, -1, Max(b,c), 1, oo
    elif b is oo and c is -oo:  # (a, oo) (-oo, d)
        br = -oo, 1, Min(a, d), -1, a, 0, oo
    elif d is oo:  # (a, b) (c, oo)
        br = -oo, -1, Min(a, c), 1, a, 0, b, -1, Max(b, c), 1, oo
    elif c is -oo:  # (a, b) (-oo, d)
        br = -oo, 1, Min(a, d), -1, a, 0, b, 1, Max(b, d), -1, oo
    elif b is oo:  # (a, oo) (c, d)
        br = -oo, -1, Min(a, c), 1, Min(a, d), -1, a, 0, oo
    elif a is -oo:  # (-oo, b) (c, d)
        br = -oo, 0, b, -1, Max(b, c), 1, Max(b, d), -1, oo
    else:  # (a, b) (c, d)
        br = (-oo, -1, Min(a, c), 1, Min(a, d), -1, a, 0, b,
            -1, Max(b, c), 1, Max(b, d), -1, oo)
    abi = []
    for _ in range(1, len(br) - 1, 2):
        a, b, i = br[_ - 1], br[_ + 1], br[_]
        if a == b:
            continue
        abi.append((a, b, i))
    if not union and not dict:
        return abi
    if union:
        rv = []
        for a, b, i in abi:
            if i == -1:
                continue
            K = (a, b)
            if not rv or a != rv[-1][1]:
                # it's a new interval
                rv.append(K)
            else:
                # it joins with the last one
                rv[-1] = rv[-1][0], K[1]
        return rv
    # dict
    rv = {}
    for a, b, i in abi:
        rv.setdefault(i, []).append((a, b))
    return rv


def _clip(A, B, k):
    """Return interval B as intervals that are covered by A (keyed
    to k) and all others keyed to -1.

    Examples
    ========

    >>> from sympy.functions.elementary.piecewise import _clip
    >>> _clip((1, 3), (2, 4), 0)
    [(2, 3, 0), (3, 4, -1)]
    """
    # pair(A, B) will give keys of -1, 0 and 1
    # we rekey A intervals to k and change all others
    # to -1
    a, b = B
    c, d = A
    A = c, d = Min(Max(c, a), b), Min(Max(d, a), b)
    p = _pair(A, B, abi=True)
    # remove infinities that weren't there at the start
    oo = S.Infinity
    if -oo not in (a, c) and p[0][0] == -oo:
        p.pop(0)
    if oo not in (b, d) and p[-1][1] == oo:
        p.pop()
    # rekey key 0 to k
    assert k != -1
    for i in range(len(p)):
        if p[i][-1] == 0:
            p[i] = p[i][:2] + (k,)
        else:
            # rekey all other keys to -1
            p[i] = p[i][:2] + (-1,)
    # combine intervals that are now at the same level
    remove = []
    for i in range(1,len(p)):
        if p[i][-1] == p[i - 1][-1]:
            remove.append(i)
            p[i - 1]=(p[i - 1][0],) + p[i][1:]
    for i in reversed(remove):
        p.pop(i)
    return [tuple(map(_mmsimp, (a, b))) + (i,) for a, b, i in p]


def _mmsimp(e):
    # simplify min/max relationships
    from sympy import simplify

    def rw(m):
        if isinstance(m, Min):
            return simplify(And(*[rw(i) for i in m.args]))
        elif isinstance(m, Max):
            return simplify(Or(*[rw(i) for i in m.args]))
        elif not m.is_Atom:
            return m.func(*[rw(a) for a in m.args])
        else:
            return m

    def wr(m):
        if isinstance(m, And):
            return Min(*[wr(i) for i in m.args])
        elif isinstance(m, Or):
            return Max(*[wr(i) for i in m.args])
        elif not m.is_Atom:
            return m.func(*[wr(a) for a in m.args])
        else:
            return m

    reps = []
    from sympy import Dummy, Number
    for m in e.atoms(Min, Max):
        nums = dict([(n, Dummy()) for n in m.atoms(Number)])
        reps.append((m, wr(rw(m.xreplace(nums))).xreplace(
            dict([(v, k) for k, v in nums.items()]))))
    return e.xreplace(dict(reps))
