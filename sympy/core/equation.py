"""
Algebraic Equations with SymPy
==============================

These tools define relations that all high school and college students would
recognize as mathematical equations. They consist of a left hand side (lhs)
and a right hand side (rhs) connected by a relation operator such as "=". At
present the "=" relation operator is the only option. The relation operator may
not be set.

This class should not be confused with the Boolean class ``Equality``
(abbreviated ``Eq``) which specifies that the equality of two objects is
``True``.

This tool applies operations to both sides of the equation simultaneously, just
as students are taught to do when attempting to isolate (solve for) a
variable. Thus the statement ``Equation/b`` yields a new equation
``Equation.lhs/b = Equation.rhs/b``

The intent is to allow using the mathematical tools in SymPy to rearrange
equations and perform algebra in a stepwise fashion. In this way more people
can successfully perform algebraic rearrangements without stumbling over
missed details such as a negative sign. This mimics the capabilities available
in [SageMath](https://www.sagemath.org/) and
[Maxima](http://maxima.sourceforge.net/).
"""


from .expr import Expr
from .basic import Basic
from .sympify import _sympify


class Equation(Basic):
    """
    This class defines an equation with a left-hand-side (lhs) and a right-
    hand-side (rhs) connected by the "=" operator (e.g. $p*V = n*R*T$).

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
    kwargs: key word arguments from this list
        - ``check=True/False`` use ``False`` to minimize computationally
           expensive calls to ``simplify`` when ``Equation`` is called
           programmatically.

    Examples
    ========
    >>> from sympy import var, Equation, Eqn, exp, log, integrate, Integral
    >>> from sympy import solve
    >>> a, b, c = var('a b c')
    >>> Equation(a,b/c)
    a = b/c
    >>> t=Eqn(a,b/c)
    >>> t
    a = b/c
    >>> t*c
    a*c = b
    >>> c*t
    a*c = b
    >>> exp(t)
    exp(a) = exp(b/c)
    >>> exp(log(t))
    a = b/c

    Utility operations
    >>> t.reversed
    b/c = a
    >>> t.swap
    b/c = a
    >>> t.lhs
    a
    >>> t.rhs
    b/c
    >>> t.as_Boolean()
    Eq(a, b/c)

    Integration can only be performed on one side at a time.
    >>> q=Eqn(a*c,b/c)
    >>> integrate(q,b,side='rhs')
    b**2/(2*c)
    >>> integrate(q,b,side='lhs')
    a*b*c

    Make a pretty statement of integration from an equation
    >>> Eqn(Integral(q.lhs,b),integrate(q,b,side='rhs'))
    Integral(a*c, b) = b**2/(2*c)

    This is duplicated by the convenience function self.integ
    >>> q.integ(b)
    Integral(a*c, b) = b**2/(2*c)

    SymPy's solvers do not understand these equations. They expect an
    expression that the solver assumes = 0. Thus to use the solver the
    equation must be rearranged so that all non-zero symbols are on one side.
    Then just the non-zero symbolic side is passed to ``solve()``.
    >>> t2 = t-t.rhs
    >>> t2
    a - b/c = 0
    >>> solve(t2.lhs,c)
    [b/a]
    """

    def __new__(cls, lhs, rhs, **kwargs):
        check = kwargs.pop('check', True)
        lhs = _sympify(lhs)
        rhs = _sympify(rhs)
        if not isinstance(lhs,Expr) or not isinstance(rhs,Expr):
            raise TypeError('lhs and rhs must be valid sympy expressions.')
        if check and lhs.is_number and rhs.is_number:
            tst=(lhs-rhs)
            if tst.is_number:
                print('Lhs & rhs are numbers. Checking statement validity.'
                      'Please be patient. This can be slow for complex '
                      'expressions...')
                tst=tst.simplify().evalf()
            if tst != 0 and tst.is_number:
                from warnings import warn
                warnstr = '\nDid you really mean to define unequal '
                warnstr += 'numbers '
                warnstr += str(lhs) + ' and ' + str(rhs)
                warnstr += ' as equal?\n'
                warnstr += 'To suppress this warning include `check=False`'
                warnstr += ' in the equation definition: '
                warnstr += '`Eqn(lhs,rhs, check=False)`.'
                warn(warnstr)
        return super().__new__(cls, lhs, rhs)

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
        from .relational import Equality
        return Equality(self.lhs, self.rhs)

    @property
    def reversed(self):
        """
        Swaps the lhs and the rhs.
        """
        return Equation(self.rhs, self.lhs, check=False)

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
        funcname = getattr(func, '__name__')
        if hasattr(expr, funcname):
            return getattr(expr, funcname)(*args, **kwargs)
        else:
            return func(expr, *args, **kwargs)

    def _applyfunc(self, func, *args, **kwargs):
        # Logic function to allow using a keyword to determine
        # which side to apply an operation/function to.

        side = kwargs.pop('Eqn_apply_side', None)
        if side == 'both':
            lhs = self._applytoexpr(self.lhs, func, *args, **kwargs)
            rhs = self._applytoexpr(self.rhs, func, *args, **kwargs)
        elif side == 'rhs':
            lhs = self.lhs
            rhs = self._applytoexpr(self.rhs, func, *args, **kwargs)
        elif side == 'lhs':
            rhs = self.rhs
            lhs = self._applytoexpr(self.lhs, func, *args, **kwargs)
        else:
            raise ValueError(
                'keyword `Eqn_apply_side` must be one of "both", "lhs" or '
                '"rhs".')
        return Equation(lhs, rhs, check=False)

    def applyfunc(self, func, *args, **kwargs):
        """
        If either side of the equation has a defined subfunction (attribute) of
        name ``func``, that will be applied instead of the global function.
        The operation is applied to both sides.
         """
        return self._applyfunc(func, *args, **kwargs, Eqn_apply_side='both')

    def applylhs(self, func, *args, **kwargs):
        """
        If lhs side of the equation has a defined subfunction (attribute) of
        name ``func``, that will be applied instead of the global function.
        The operation is applied to only the lhs.
        """
        return self._applyfunc(func, *args, **kwargs, Eqn_apply_side='lhs')

    def applyrhs(self, func, *args, **kwargs):
        """
        If rhs side of the equation has a defined subfunction (attribute) of
        name ``func``, that will be applied instead of the global function.
        The operation is applied to only the rhs.
        """
        return self._applyfunc(func, *args, **kwargs, Eqn_apply_side='rhs')

    #####
    # Overrides of binary math operations
    #####

    @classmethod
    def _binary_op(cls, a, b, opfunc_ab):
        if isinstance(a, Equation) and not isinstance(b, Equation):
            return Equation(opfunc_ab(a.lhs, b), opfunc_ab(a.rhs, b),
                            check=False)
        elif isinstance(b, Equation) and not isinstance(a, Equation):
            return Equation(opfunc_ab(a, b.lhs), opfunc_ab(a, b.rhs),
                            check=False)
        elif isinstance(a, Equation) and isinstance(b, Equation):
            return Equation(opfunc_ab(a.lhs, b.lhs), opfunc_ab(a.rhs, b.rhs),
                            check=False)
        else:
            raise TypeError('One of a or b should be an equation')

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
            *args, **kwargs), check=False)

    def simplify(self, *args, **kwargs):
        return self._eval_simplify(*args, **kwargs)

    def _eval_simplify(self, *args, **kwargs):
        return Equation(self.lhs.simplify(*args, **kwargs), self.rhs.simplify(
            *args, **kwargs), check=False)

    def _eval_factor(self, *args, **kwargs):
        # TODO: cancel out factors common to both sides.
        return Equation(self.lhs.factor(*args, **kwargs), self.rhs.factor(
            *args, **kwargs), check=False)

    def factor(self, *args, **kwargs):
        return self._eval_factor(*args, **kwargs)

    def _eval_collect(self, *args, **kwargs):
        from sympy.simplify.radsimp import collect
        return Equation(collect(self.lhs, *args, **kwargs),
                        collect(self.rhs, *args, **kwargs), check=False)

    def collect(self, *args, **kwargs):
        return self._eval_collect(*args, **kwargs)

    def evalf(self, *args, **kwargs):
        return Equation(self.lhs.evalf(*args, **kwargs),
                        self.rhs.evalf(*args, **kwargs), check=False)

    def _eval_derivative(self, *args, **kwargs):
        # TODO Find why diff and Derivative do not appear to pass through
        #  kwargs to this. Since we cannot set evaluation of lhs manually
        #  try to be intelligent about when to do it.
        from sympy.core.function import Derivative
        eval_lhs = False
        if not (isinstance(self.lhs, Derivative)):
            for sym in args:
                if sym in self.lhs.free_symbols and not (
                _sympify(sym).is_number):
                    eval_lhs = True
        return Equation(self.lhs.diff(*args, **kwargs, evaluate=eval_lhs),
                        self.rhs.diff(*args, **kwargs), check=False)

    def integ(self, *args, **kwargs):
        """
        This function is a convenience function that returns a new equation
        consisting of an unevaluated integral of the lhs as the new lhs and
        the result of integrating the rhs as the new rhs.
        """
        from sympy.integrals.integrals import Integral
        return Equation(Integral(self.lhs, *args, **kwargs),
                        self.rhs.integrate(*args, **kwargs), check=False)

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


Eqn = Equation
