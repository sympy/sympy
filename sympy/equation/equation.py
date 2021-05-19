from sympy.core.expr import Expr
from sympy.core.basic import Basic
from sympy.core.evalf import EvalfMixin
from sympy.core.sympify import _sympify
import functools


class Equation(Basic, EvalfMixin):
    """
    This class defines an equation with a left-hand-side (lhs) and a right-
    hand-side (rhs) connected by the "=" operator (e.g. p*V = n*R*T). In
    interactive environments (e.g. Jupyter) with pretty printing available
    this displays in standard mathematical notation.

    Explanation
    ===========
    This class defines relations that all high school and college students
    would recognize as mathematical equations. Only the "=" relation
    operator is recognized.

    This class is intended to allow using the mathematical tools in SymPy to
    rearrange equations and perform algebra in a stepwise fashion. In this
    way more people can successfully perform algebraic rearrangements without
    stumbling over missed details such as a negative sign.

    Create an equation with the call ``Equation(lhs,rhs)``, where ``lhs`` and
    ``rhs`` are any valid Sympy expression. ``Eqn(...)`` is a synonym for
    ``Equation(...)``.

    In a text only environment or if using ``print(Eqn(...))`` equations
    will be printed as ``Equation(lhs,rhs)``. Thus all the docstring examples
    below do not display as ``lhs = rhs``, but as ``Equation(lhs,rhs)``.

    Parameters
    ==========
    lhs: sympy expression, ``class Expr``.
    rhs: sympy expression, ``class Expr``.
    kwargs:

    Examples
    ========
    >>> from sympy import var, Equation, Eqn, exp, log, integrate, Integral
    >>> from sympy import factor, collect
    >>> from sympy import solve
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

    Differentiation should only be applied to one side at a time to avoid
    unexpected results (e.g. 0 = non-zero expression).
    >>> q=Eqn(a*c, b/c**2)
    >>> q
    Equation(a*c, b/c**2)
    >>> Eqn(q.lhs.diff(c),q.rhs.diff(c))
    Equation(a, -2*b/c**3)
    >>> q.dorhs.diff(c).dolhs.diff(c)
    Equation(a, -2*b/c**3)

    Integration can only be performed on one side at a time.
    Make a pretty statement of integration from an equation
    >>> Eqn(Integral(q.lhs,b),integrate(q.rhs,b))
    Equation(Integral(a*c, b), b**2/(2*c**2))

    Integration of each side with respect to different variables
    >>> q.dorhs.integrate(b).dolhs.integrate(a)
    Equation(a**2*c/2, b**2/(2*c**2))

    SymPy's solvers do not understand these equations. They expect an
    expression that the solver assumes = 0. Thus to use the solver the
    equation must be rearranged so that all non-zero symbols are on one side.
    Then just the non-zero symbolic side is passed to ``solve()``.
    >>> t2 = t-t.rhs
    >>> t2
    Equation(a - b/c, 0)
    >>> solve(t2.lhs,c)
    [b/a]
    """

    def __new__(cls, lhs, rhs, **kwargs):
        lhs = _sympify(lhs)
        rhs = _sympify(rhs)
        if not isinstance(lhs, Expr) or not isinstance(rhs, Expr):
            raise TypeError('lhs and rhs must be valid sympy expressions.')
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

    def as_Boolean(self, **kwargs):
        """
        Converts the equation to an `Equality`.
        """
        from ..core.relational import Equality
        return Equality(self.lhs, self.rhs, **kwargs)

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
        from ..core.relational import Equality
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
            func = None
            if self.side in ('rhs', 'both'):
                func = getattr(self.eqn.rhs, name, None)
            else:
                func = getattr(self.eqn.lhs, name, None)
            if func is None:
                raise AttributeError('Expressions in the equation have no '
                                     'attribute `'+str(name)+'`. Try `.apply('
                                     + str(name)+', *args)` or '
                                     'pass the equation as a parameter to `'
                                     + str(name)+'()`.')
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


Eqn = Equation
