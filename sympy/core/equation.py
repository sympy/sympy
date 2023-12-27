from sympy.core.basic import Basic
from sympy.core.evalf import EvalfMixin

class Equation(Basic, EvalfMixin):
    """
    This class defines an equation with a left-hand-side (lhs) and a right-
    hand-side (rhs) connected by the "=" operator (e.g. `p*V = n*R*T`).

    Explanation
    ===========
    This class defines relations that all high school and college students
    would recognize as mathematical equations. At present only the "=" relation
    operator is recognized.

    This class is intended to allow using the mathematical tools in SymPy to
    rearrange equations and perform algebra in a stepwise fashion. In this
    way more people can successfully perform algebraic rearrangements without
    stumbling over missed details such as a negative sign.

    Create an equation with the call ``Equation(lhs,rhs)``, where ``lhs`` and
    ``rhs`` are any valid Sympy expression. ``Eqn(...)`` is a synonym for
    ``Equation(...)``.

    Parameters
    ==========
    lhs: sympy expression, ``class Expr``.
    rhs: sympy expression, ``class Expr``.
    kwargs:

    Examples
    ========
    NOTE: If used with `algebra_with_sympy`
    (https://github.com/gutow/Algebra_with_Sympy) you can get human-readable
    output.
    >>> from sympy import var, Equation, Eqn, exp, log, diff
    >>> from sympy import integrate, Integral
    >>> a, b, c, x = var('a b c x')
    >>> Equation(a,b/c)
    Equation(a, b/c)
    >>> t=Eqn(a,b/c)
    >>> t
    Equation(a, b/c)
    >>> t*c
    Equation(a*c, b)
    >>> c*t
    Equation(a*c, b)
    >>> exp(t)
    Equation(exp(a), exp(b/c))
    >>> exp(log(t))
    Equation(a, b/c)

    Simplification and Expansion
    >>> from sympy import simplify, expand, factor, collect
    >>> f = Eqn(x**2 - 1, c)
    >>> f
    Equation(x**2 - 1, c)
    >>> f/(x+1)
    Equation((x**2 - 1)/(x + 1), c/(x + 1))
    >>> (f/(x+1)).simplify()
    Equation(x - 1, c/(x + 1))
    >>> simplify(f/(x+1))
    Equation(x - 1, c/(x + 1))
    >>> (f/(x+1)).expand()
    Equation(x**2/(x + 1) - 1/(x + 1), c/(x + 1))
    >>> expand(f/(x+1))
    Equation(x**2/(x + 1) - 1/(x + 1), c/(x + 1))
    >>> factor(f)
    Equation((x - 1)*(x + 1), c)
    >>> f.factor()
    Equation((x - 1)*(x + 1), c)
    >>> f2 = f+a*x**2+b*x +c
    >>> f2
    Equation(a*x**2 + b*x + c + x**2 - 1, a*x**2 + b*x + 2*c)
    >>> collect(f2,x)
    Equation(b*x + c + x**2*(a + 1) - 1, a*x**2 + b*x + 2*c)

    Apply operation to only one side
    >>> poly = Eqn(a*x**2 + b*x + c*x**2, a*x**3 + b*x**3 + c*x)
    >>> poly.applyrhs(factor,x)
    Equation(a*x**2 + b*x + c*x**2, x*(c + x**2*(a + b)))
    >>> poly.applylhs(factor)
    Equation(x*(a*x + b + c*x), a*x**3 + b*x**3 + c*x)
    >>> poly.applylhs(collect,x)
    Equation(b*x + x**2*(a + c), a*x**3 + b*x**3 + c*x)

    ``.apply...`` also works with user defined python functions
    >>> def addsquare(eqn):
    ...     return eqn+eqn**2
    ...
    >>> t.apply(addsquare)
    Equation(a**2 + a, b**2/c**2 + b/c)
    >>> t.applyrhs(addsquare)
    Equation(a, b**2/c**2 + b/c)
    >>> t.apply(addsquare, side = 'rhs')
    Equation(a, b**2/c**2 + b/c)
    >>> t.applylhs(addsquare)
    Equation(a**2 + a, b/c)
    >>> addsquare(t)
    Equation(a**2 + a, b**2/c**2 + b/c)

    Inaddition to ``.apply...`` there is also the less general ``.do``,
    ``.dolhs``, ``.dorhs``, which only works for operations defined on the
    ``Expr`` class (e.g.``.collect(), .factor(), .expand()``, etc...).
    >>> poly.dolhs.collect(x)
    Equation(b*x + x**2*(a + c), a*x**3 + b*x**3 + c*x)
    >>> poly.dorhs.collect(x)
    Equation(a*x**2 + b*x + c*x**2, c*x + x**3*(a + b))
    >>> poly.do.collect(x)
    Equation(b*x + x**2*(a + c), c*x + x**3*(a + b))
    >>> poly.dorhs.factor()
    Equation(a*x**2 + b*x + c*x**2, x*(a*x**2 + b*x**2 + c))

    ``poly.do.exp()`` or other sympy math functions will raise an error.

    Rearranging an equation (simple example made complicated as illustration)
    >>> p, V, n, R, T = var('p V n R T')
    >>> eq1=Eqn(p*V,n*R*T)
    >>> eq1
    Equation(V*p, R*T*n)
    >>> eq2 =eq1/V
    >>> eq2
    Equation(p, R*T*n/V)
    >>> eq3 = eq2/R/T
    >>> eq3
    Equation(p/(R*T), n/V)
    >>> eq4 = eq3*R/p
    >>> eq4
    Equation(1/T, R*n/(V*p))
    >>> 1/eq4
    Equation(T, V*p/(R*n))
    >>> eq5 = 1/eq4 - T
    >>> eq5
    Equation(0, -T + V*p/(R*n))

    Substitution (#'s and units)
    >>> L, atm, mol, K = var('L atm mol K', positive=True, real=True) # units
    >>> eq2.subs({R:0.08206*L*atm/mol/K,T:273*K,n:1.00*mol,V:24.0*L})
    Equation(p, 0.9334325*atm)
    >>> eq2.subs({R:0.08206*L*atm/mol/K,T:273*K,n:1.00*mol,V:24.0*L}).evalf(4)
    Equation(p, 0.9334*atm)

    Substituting an equation into another equation:
    >>> P, P1, P2, A1, A2, E1, E2 = var("P, P1, P2, A1, A2, E1, E2")
    >>> eq1 = Eqn(P, P1 + P2)
    >>> eq2 = Eqn(P1 / (A1 * E1), P2 / (A2 * E2))
    >>> P1_val = (eq1 - P2).swap
    >>> P1_val
    Equation(P1, P - P2)
    >>> eq2 = eq2.subs(P1_val)
    >>> eq2
    Equation((P - P2)/(A1*E1), P2/(A2*E2))

    Combining equations (Math with equations: lhs with lhs and rhs with rhs)
    >>> q = Eqn(a*c, b/c**2)
    >>> q
    Equation(a*c, b/c**2)
    >>> t
    Equation(a, b/c)
    >>> q+t
    Equation(a*c + a, b/c + b/c**2)
    >>> q/t
    Equation(c, 1/c)
    >>> t**q
    Equation(a**(a*c), (b/c)**(b/c**2))

    Utility operations
    >>> t.reversed
    Equation(b/c, a)
    >>> t.swap
    Equation(b/c, a)
    >>> t.lhs
    a
    >>> t.rhs
    b/c
    >>> t.as_Boolean()
    Eq(a, b/c)

    `.check()` convenience method for `.as_Boolean().simplify()`
    >>> from sympy import I, pi
    >>> Equation(pi*(I+2), pi*I+2*pi).check()
    True
    >>> Eqn(a,a+1).check()
    False

    Differentiation
    Differentiation is applied to both sides if the wrt variable appears on
    both sides.
    >>> q=Eqn(a*c, b/c**2)
    >>> q
    Equation(a*c, b/c**2)
    >>> diff(q,b)
    Equation(Derivative(a*c, b), c**(-2))
    >>> diff(q,c)
    Equation(a, -2*b/c**3)
    >>> diff(log(q),b)
    Equation(Derivative(log(a*c), b), 1/b)
    >>> diff(q,c,2)
    Equation(Derivative(a, c), 6*b/c**4)

    If you specify multiple differentiation all at once the assumption
    is order of differentiation matters and the lhs will not be
    evaluated.
    >>> diff(q,c,b)
    Equation(Derivative(a*c, b, c), -2/c**3)

    To overcome this specify the order of operations.
    >>> diff(diff(q,c),b)
    Equation(Derivative(a, b), -2/c**3)

    But the reverse order returns an unevaulated lhs (a may depend on b).
    >>> diff(diff(q,b),c)
    Equation(Derivative(a*c, b, c), -2/c**3)

    Integration can only be performed on one side at a time.
    >>> q=Eqn(a*c,b/c)
    >>> integrate(q,b,side='rhs')
    b**2/(2*c)
    >>> integrate(q,b,side='lhs')
    a*b*c

    Make a pretty statement of integration from an equation
    >>> Eqn(Integral(q.lhs,b),integrate(q,b,side='rhs'))
    Equation(Integral(a*c, b), b**2/(2*c))

    Integration of each side with respect to different variables
    >>> q.dorhs.integrate(b).dolhs.integrate(a)
    Equation(a**2*c/2, b**2/(2*c))
    """

    def __new__(cls, lhs, rhs, **kwargs):
        from sympy.core.sympify import _sympify
        from sympy.core.expr import Expr
        lhs = _sympify(lhs)
        rhs = _sympify(rhs)
        if not isinstance(lhs, Expr) or not isinstance(rhs, Expr):
            raise TypeError('lhs and rhs must be valid sympy expressions.')
        return super().__new__(cls, lhs, rhs)

    def _get_eqn_name(self):
        """
        Tries to find the python string name that refers to the equation. In
        IPython environments (IPython, Jupyter, etc...) looks in the user_ns.
        If not in an IPython environment looks in __main__.
        :return: string value if found or empty string.
        """
        import __main__ as shell
        for k in dir(shell):
            item = getattr(shell, k)
            if isinstance(item, Equation):
                if item == self and not k.startswith('_'):
                    return k
        return ''

    @property
    def lhs(self):
        """
        Returns the lhs of the equation.
        """
        return self.args[0]

    @property
    def rhs(self):
        """
        Returns the rhs of the equation.
        """
        return self.args[1]

    def as_Boolean(self):
        """
        Converts the equation to an Equality.
        """
        from sympy import Equality
        return Equality(self.lhs, self.rhs)

    def check(self, **kwargs):
        """
        Forces simplification and casts as `Equality` to check validity.
        Parameters
        ----------
        kwargs any appropriate for `Equality`.

        Returns
        -------
        True, False or an unevaluated `Equality` if truth cannot be determined.
        """
        from sympy.core.relational import Equality
        return Equality(self.lhs, self.rhs, **kwargs).simplify()

    @property
    def reversed(self):
        """
        Swaps the lhs and the rhs.
        """
        return Equation(self.rhs, self.lhs)

    @property
    def swap(self):
        """
        Synonym for `.reversed`
        """
        return self.reversed

    def _applytoexpr(self, expr, func, *args, **kwargs):
        # Applies a function to an expression checking whether there
        # is a specialized version associated with the particular type of
        # expression. Errors will be raised if the function cannot be
        # applied to an expression.
        funcname = getattr(func, '__name__', None)
        if funcname is not None:
            localfunc = getattr(expr, funcname, None)
            if localfunc is not None:
                return localfunc(*args, **kwargs)
        return func(expr, *args, **kwargs)

    def apply(self, func, *args, side='both', **kwargs):
        """
        Apply an operation/function/method to the equation returning the
        resulting equation.

        Parameters
        ==========

        func: object
            object to apply usually a function

        args: as necessary for the function

        side: 'both', 'lhs', 'rhs', optional
            Specifies which side of the equation the operation will be applied
            to. Default is 'both'.

        kwargs: as necessary for the function
         """
        lhs = self.lhs
        rhs = self.rhs
        if side in ('both', 'lhs'):
            lhs = self._applytoexpr(self.lhs, func, *args, **kwargs)
        if side in ('both', 'rhs'):
            rhs = self._applytoexpr(self.rhs, func, *args, **kwargs)
        return Equation(lhs, rhs)

    def applylhs(self, func, *args, **kwargs):
        """
        If lhs side of the equation has a defined subfunction (attribute) of
        name ``func``, that will be applied instead of the global function.
        The operation is applied to only the lhs.
        """
        return self.apply(func, *args, **kwargs, side='lhs')

    def applyrhs(self, func, *args, **kwargs):
        """
        If rhs side of the equation has a defined subfunction (attribute) of
        name ``func``, that will be applied instead of the global function.
        The operation is applied to only the rhs.
        """
        return self.apply(func, *args, **kwargs, side='rhs')

    class _sides:
        """
        Helper class for the `.do.`, `.dolhs.`, `.dorhs.` syntax for applying
        submethods of expressions.
        """

        def __init__(self, eqn, side='both'):
            self.eqn = eqn
            self.side = side

        def __getattr__(self, name):
            import functools
            func = None
            if self.side in ('rhs', 'both'):
                func = getattr(self.eqn.rhs, name, None)
            else:
                func = getattr(self.eqn.lhs, name, None)
            if func is None:
                raise AttributeError('Expressions in the equation have no '
                                     'attribute `' + str(
                    name) + '`. Try `.apply('
                                     + str(name) + ', *args)` or '
                                                   'pass the equation as a parameter to `'
                                     + str(name) + '()`.')
            return functools.partial(self.eqn.apply, func, side=self.side)

    @property
    def do(self):
        return self._sides(self, side='both')

    @property
    def dolhs(self):
        return self._sides(self, side='lhs')

    @property
    def dorhs(self):
        return self._sides(self, side='rhs')

    def _eval_rewrite(self, rule, args, **kwargs):
        """Return Equation(L, R) as Equation(L - R, 0) or as L - R.

        Parameters
        ==========

        evaluate : bool, optional
            Control the evaluation of the result. If `evaluate=None` then
            terms in L and R will not cancel but they will be listed in
            canonical order; otherwise non-canonical args will be returned.
            Default to True.

        eqn : bool, optional
            Control the returned type. If `eqn=True`, then Equation(L - R, 0)
            is returned. Otherwise, the L - R symbolic expression is returned.
            Default to True.

        Examples
        ========
        >>> from sympy import Add
        >>> from sympy.abc import b, x
        >>> from sympy import Equation
        >>> eq = Equation(x + b, x - b)
        >>> eq.rewrite(Add)
        Equation(2*b, 0)
        >>> eq.rewrite(Add, evaluate=None).lhs.args
        (b, b, x, -x)
        >>> eq.rewrite(Add, evaluate=False).lhs.args
        (b, x, b, -x)
        >>> eq.rewrite(Add, eqn=False)
        2*b
        >>> eq.rewrite(Add, eqn=False, evaluate=False).args
        (b, x, b, -x)
        """
        from sympy import Add
        from sympy.core.add import _unevaluated_Add
        if rule == Add:
            # NOTE: the code about `evaluate` is very similar to
            # sympy.core.relational.Equality._eval_rewrite_as_Add
            eqn = kwargs.pop("eqn", True)
            evaluate = kwargs.get('evaluate', True)
            L, R = args
            if evaluate:
                # allow cancellation of args
                expr = L - R
            else:
                args = Add.make_args(L) + Add.make_args(-R)
                if evaluate is None:
                    # no cancellation, but canonical
                    expr = _unevaluated_Add(*args)
                else:
                    # no cancellation, not canonical
                    expr = Add._from_args(args)
            if eqn:
                return self.func(expr, 0)
            return expr

    def subs(self, *args, **kwargs):
        """Substitutes old for new in an equation after sympifying args.

        `args` is either:

        * one or more arguments of type `Equation(old, new)`.
        * two arguments, e.g. foo.subs(old, new)
        * one iterable argument, e.g. foo.subs(iterable). The iterable may be:

            - an iterable container with (old, new) pairs. In this case the
              replacements are processed in the order given with successive
              patterns possibly affecting replacements already made.
            - a dict or set whose key/value items correspond to old/new pairs.
              In this case the old/new pairs will be sorted by op count and in
              case of a tie, by number of args and the default_sort_key. The
              resulting sorted list is then processed as an iterable container
              (see previous).

        If the keyword ``simultaneous`` is True, the subexpressions will not be
        evaluated until all the substitutions have been made.

        Please, read ``help(Expr.subs)`` for more examples.

        Examples
        ========

        >>> from sympy.abc import a, b, c, x
        >>> from sympy import Equation
        >>> eq = Equation(x + a, b * c)

        Substitute a single value:

        >>> eq.subs(b, 4)
        Equation(a + x, 4*c)

        Substitute a multiple values:

        >>> eq.subs([(a, 2), (b, 4)])
        Equation(x + 2, 4*c)
        >>> eq.subs({a: 2, b: 4})
        Equation(x + 2, 4*c)

        Substitute an equation into another equation:

        >>> eq2 = Equation(x + a, 4)
        >>> eq.subs(eq2)
        Equation(4, b*c)

        Substitute multiple equations into another equation:

        >>> eq1 = Equation(x + a + b + c, x * a * b * c)
        >>> eq2 = Equation(x + a, 4)
        >>> eq3 = Equation(b, 5)
        >>> eq1.subs(eq2, eq3)
        Equation(c + 9, 5*a*c*x)

        """
        new_args = args
        if all(isinstance(a, self.func) for a in args):
            new_args = [{a.args[0]: a.args[1] for a in args}]
        elif (len(args) == 1) and all(isinstance(a, self.func) for a in
                                      args[0]):
            raise TypeError("You passed into `subs` a list of elements of "
                            "type `Equation`, but this is not supported. Please, consider "
                            "unpacking the list with `.subs(*eq_list)` or select your "
                            "equations from the list and use `.subs(eq_list[0], eq_list["
                            "2], ...)`.")
        elif any(isinstance(a, self.func) for a in args):
            raise ValueError("`args` contains one or more Equation and some "
                             "other data type. This mode of operation is not supported. "
                             "Please, read `subs` documentation to understand how to "
                             "use it.")
        return super().subs(*new_args, **kwargs)

    #####
    # Overrides of binary math operations
    #####

    @classmethod
    def _binary_op(cls, a, b, opfunc_ab):
        if isinstance(a, Equation) and not isinstance(b, Equation):
            return Equation(opfunc_ab(a.lhs, b), opfunc_ab(a.rhs, b))
        elif isinstance(b, Equation) and not isinstance(a, Equation):
            return Equation(opfunc_ab(a, b.lhs), opfunc_ab(a, b.rhs))
        elif isinstance(a, Equation) and isinstance(b, Equation):
            return Equation(opfunc_ab(a.lhs, b.lhs), opfunc_ab(a.rhs, b.rhs))
        else:
            return NotImplemented

    def __add__(self, other):
        return self._binary_op(self, other, lambda a, b: a + b)

    def __radd__(self, other):
        return self._binary_op(other, self, lambda a, b: a + b)

    def __mul__(self, other):
        return self._binary_op(self, other, lambda a, b: a * b)

    def __rmul__(self, other):
        return self._binary_op(other, self, lambda a, b: a * b)

    def __sub__(self, other):
        return self._binary_op(self, other, lambda a, b: a - b)

    def __rsub__(self, other):
        return self._binary_op(other, self, lambda a, b: a - b)

    def __truediv__(self, other):
        return self._binary_op(self, other, lambda a, b: a / b)

    def __rtruediv__(self, other):
        return self._binary_op(other, self, lambda a, b: a / b)

    def __mod__(self, other):
        return self._binary_op(self, other, lambda a, b: a % b)

    def __rmod__(self, other):
        return self._binary_op(other, self, lambda a, b: a % b)

    def __pow__(self, other):
        return self._binary_op(self, other, lambda a, b: a ** b)

    def __rpow__(self, other):
        return self._binary_op(other, self, lambda a, b: a ** b)

    def _eval_power(self, other):
        return self.__pow__(other)

    #####
    # Operation helper functions
    #####
    def expand(self, *args, **kwargs):
        return Equation(self.lhs.expand(*args, **kwargs), self.rhs.expand(
            *args, **kwargs))

    def simplify(self, *args, **kwargs):
        return self._eval_simplify(*args, **kwargs)

    def _eval_simplify(self, *args, **kwargs):
        return Equation(self.lhs.simplify(*args, **kwargs), self.rhs.simplify(
            *args, **kwargs))

    def _eval_factor(self, *args, **kwargs):
        # TODO: cancel out factors common to both sides.
        return Equation(self.lhs.factor(*args, **kwargs), self.rhs.factor(
            *args, **kwargs))

    def factor(self, *args, **kwargs):
        return self._eval_factor(*args, **kwargs)

    def _eval_collect(self, *args, **kwargs):
        from sympy.simplify.radsimp import collect
        return Equation(collect(self.lhs, *args, **kwargs),
                        collect(self.rhs, *args, **kwargs))

    def collect(self, *args, **kwargs):
        return self._eval_collect(*args, **kwargs)

    def evalf(self, *args, **kwargs):
        return Equation(self.lhs.evalf(*args, **kwargs),
                        self.rhs.evalf(*args, **kwargs))

    n = evalf

    def _eval_derivative(self, *args, **kwargs):
        # TODO Find why diff and Derivative do not appear to pass through
        #  kwargs to this. Since we cannot set evaluation of lhs manually
        #  try to be intelligent about when to do it.
        from sympy.core.function import Derivative
        from sympy.core.sympify import _sympify
        eval_lhs = False
        if not (isinstance(self.lhs, Derivative)):
            for sym in args:
                if sym in self.lhs.free_symbols and not (
                    _sympify(sym).is_number):
                    eval_lhs = True
        return Equation(self.lhs.diff(*args, **kwargs, evaluate=eval_lhs),
                        self.rhs.diff(*args, **kwargs))

    def _eval_Integral(self, *args, **kwargs):
        side = kwargs.pop('side', None)  # Could not seem to pass values for
        # `evaluate` through to here.
        if side is None:
            raise ValueError('You must specify `side="lhs"` or `side="rhs"` '
                             'when integrating an Equation')
        else:
            try:
                return (getattr(self, side).integrate(*args, **kwargs))
            except AttributeError:
                raise AttributeError('`side` must equal "lhs" or "rhs".')

    #####
    # Output helper functions
    #####
    def __repr__(self):
        repstr = 'Equation(%s, %s)' % (
        self.lhs.__repr__(), self.rhs.__repr__())
        # if algwsym_config.output.human_text:
        #     return self.__str__()
        return repstr

    def _latex(self, printer):
        tempstr = ''
        tempstr += printer._print(self.lhs)
        tempstr += '='
        tempstr += printer._print(self.rhs)
        return tempstr

    def __str__(self):
        tempstr = ''
        tempstr += str(self.lhs) + ' = ' + str(self.rhs)
        return tempstr


Eqn = Equation
