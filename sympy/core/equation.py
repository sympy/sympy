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
from .evalf import EvalfMixin
from .sympify import _sympify
import functools


class Equation(Basic, EvalfMixin):
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
        - ``check=True/False``. Defaults to ``False`` to minimize
            computationally expensive calls to ``simplify`` when
            ``Equation`` is called programmatically. ``check=True``
            produces a warning if the lhs and rhs evaluate to numbers and
            are not equal.

    Examples
    ========
    >>> from sympy import var, Equation, Eqn, exp, log, integrate, Integral
    >>> from sympy import simplify, collect, expand, factor, diff
    >>> from sympy import solve
    >>> a, b, c, x = var('a b c x')
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

    Simplification and Expansion
    >>> f = Eqn(x**2 - 1, c)
    >>> f
    x**2 - 1 = c
    >>> f/(x+1)
    (x**2 - 1)/(x + 1) = c/(x + 1)
    >>> (f/(x+1)).simplify()
    x - 1 = c/(x + 1)
    >>> simplify(f/(x+1))
    x - 1 = c/(x + 1)
    >>> (f/(x+1)).expand()
    x**2/(x + 1) - 1/(x + 1) = c/(x + 1)
    >>> expand(f/(x+1))
    x**2/(x + 1) - 1/(x + 1) = c/(x + 1)
    >>> factor(f)
    (x - 1)*(x + 1) = c
    >>> f.factor()
    (x - 1)*(x + 1) = c
    >>> f2 = f+a*x**2+b*x +c
    >>> f2
    a*x**2 + b*x + c + x**2 - 1 = a*x**2 + b*x + 2*c
    >>> collect(f2,x)
    b*x + c + x**2*(a + 1) - 1 = a*x**2 + b*x + 2*c

    Apply operation to only one side
    >>> poly = Eqn(a*x**2 + b*x + c*x**2, a*x**3 + b*x**3 + c*x)
    >>> poly.applyrhs(factor,x)
    a*x**2 + b*x + c*x**2 = x*(c + x**2*(a + b))
    >>> poly.applylhs(factor)
    x*(a*x + b + c*x) = a*x**3 + b*x**3 + c*x
    >>> poly.applylhs(collect,x)
    b*x + x**2*(a + c) = a*x**3 + b*x**3 + c*x

    ``.apply...`` also works with user defined python functions
    >>> def addsquare(expr):
    ...     return expr+expr**2
    ...
    >>> t.apply(addsquare)
    a**2 + a = b**2/c**2 + b/c
    >>> t.applyrhs(addsquare)
    a = b**2/c**2 + b/c
    >>> t.apply(addsquare, side = 'rhs')
    a = b**2/c**2 + b/c
    >>> t.applylhs(addsquare)
    a**2 + a = b/c
    >>> addsquare(t)
    a**2 + a = b**2/c**2 + b/c

    Inaddition to ``.apply...`` there is also the less general ``.do``,
    ``.dolhs``, ``.dorhs``, which only works for operations defined on the
    ``Expr`` class (e.g.``.collect(), .factor(), .expand()``, etc...).
    >>> poly.dolhs.collect(x)
    b*x + x**2*(a + c) = a*x**3 + b*x**3 + c*x
    >>> poly.dorhs.collect(x)
    a*x**2 + b*x + c*x**2 = c*x + x**3*(a + b)
    >>> poly.do.collect(x)
    b*x + x**2*(a + c) = c*x + x**3*(a + b)
    >>> poly.dorhs.factor()
    a*x**2 + b*x + c*x**2 = x*(a*x**2 + b*x**2 + c)

    ``poly.do.exp()`` or other sympy math functions will raise an error.

    Rearranging an equation (simple example made complicated as illustration)
    >>> p, V, n, R, T = var('p V n R T')
    >>> eq1=Eqn(p*V,n*R*T)
    >>> eq1
    V*p = R*T*n
    >>> eq2 =eq1/V
    >>> eq2
    p = R*T*n/V
    >>> eq3 = eq2/R/T
    >>> eq3
    p/(R*T) = n/V
    >>> eq4 = eq3*R/p
    >>> eq4
    1/T = R*n/(V*p)
    >>> 1/eq4
    T = V*p/(R*n)
    >>> eq5 = 1/eq4 - T
    >>> eq5
    0 = -T + V*p/(R*n)

    Substitution (#'s and units)
    >>> L, atm, mol, K = var('L atm mol K', positive=True, real=True) # units
    >>> eq2.subs({R:0.08206*L*atm/mol/K,T:273*K,n:1.00*mol,V:24.0*L})
    p = 0.9334325*atm
    >>> eq2.subs({R:0.08206*L*atm/mol/K,T:273*K,n:1.00*mol,V:24.0*L}).evalf(4)
    p = 0.9334*atm

    Combining equations (Math with equations: lhs with lhs and rhs with rhs)
    >>> q = Eqn(a*c, b/c**2)
    >>> q
    a*c = b/c**2
    >>> t
    a = b/c
    >>> q+t
    a*c + a = b/c + b/c**2
    >>> q/t
    c = 1/c
    >>> t**q
    a**(a*c) = (b/c)**(b/c**2)

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

    Differentiation
    Differentiation is applied to both sides if the wrt variable appears on
    both sides.
    >>> q=Eqn(a*c, b/c**2)
    >>> q
    a*c = b/c**2
    >>> diff(q,b)
    Derivative(a*c, b) = c**(-2)
    >>> diff(q,c)
    a = -2*b/c**3
    >>> diff(log(q),b)
    Derivative(log(a*c), b) = 1/b
    >>> diff(q,c,2)
    Derivative(a, c) = 6*b/c**4

    If you specify multiple differentiation all at once the assumption
    is order of differentiation matters and the lhs will not be
    evaluated.
    >>> diff(q,c,b)
    Derivative(a*c, b, c) = -2/c**3

    To overcome this specify the order of operations.
    >>> diff(diff(q,c),b)
    Derivative(a, b) = -2/c**3

    But the reverse order returns an unevaulated lhs (a may depend on b).
    >>> diff(diff(q,b),c)
    Derivative(a*c, b, c) = -2/c**3

    Integration can only be performed on one side at a time.
    >>> q=Eqn(a*c,b/c)
    >>> integrate(q,b,side='rhs')
    b**2/(2*c)
    >>> integrate(q,b,side='lhs')
    a*b*c

    Make a pretty statement of integration from an equation
    >>> Eqn(Integral(q.lhs,b),integrate(q,b,side='rhs'))
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
        check = kwargs.pop('check', False)
        lhs = _sympify(lhs)
        rhs = _sympify(rhs)
        if not isinstance(lhs,Expr) or not isinstance(rhs,Expr):
            raise TypeError('lhs and rhs must be valid sympy expressions.')
        if check:
            tst=(lhs-rhs)
            tst=tst.simplify().evalf()
            if tst != 0 and tst.is_number:
                raise ValueError('`lhs-rhs` evaluates to a non-zero number.')
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

        side: 'both', 'lhs', 'rhs', optional
            Specifies which side of the equation the operation will be applied
            to. Default is 'both'.

         """
        lhs = self.lhs
        rhs = self.rhs
        if side in ('both', 'lhs'):
            lhs = self._applytoexpr(self.lhs, func, *args, **kwargs)
        if side in ('both', 'rhs'):
            rhs = self._applytoexpr(self.rhs, func, *args, **kwargs)
        return Equation(lhs, rhs, check=False)

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
        def __init__(self,eqn, side='both'):
            self.eqn = eqn
            self.side = side

        def __getattr__(self, name):
            func = None
            if self.side in ('rhs', 'both'):
                func = getattr(self.eqn.rhs, name, None)
            else:
                func = getattr(self.eqn.lhs, name, None)
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

    n = evalf

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
