from sympy.core import Expr, S, Tuple
from sympy.core.sympify import _sympify
from sympy.core.symbol import Str
from sympy.sets import ProductSet

__all__ = ['Map', 'InverseMap', 'IdentityMap', 'AppliedMap']

class Map(Expr):
    """
    A class for general mathematical mappings

    Parameters
    ==========

    parameters : tuple
        Parameters of the map

    domain, codomain : Set

    name : str
        Name of the map

    Examples
    ========

    >>> from sympy import symbols
    >>> from sympy.map import Map
    >>> x, y = symbols('x y')

    >>> class F(Map):
    ...     def eval(self, x):
    ...         a, = self.parameters
    ...         return a*x

    >>> f1 = F(parameters=(x,), name='f')
    >>> f1.name
    'f'
    >>> f1.parameters
    (x,)
    >>> f1(y).doit()
    x*y

    >>> f2 = f1.subs(x, 2)
    >>> f2.parameters
    (2,)
    >>> f2(y).doit()
    2*y

    """
    def __new__(
        cls, parameters=(), domain=S.UniversalSet, codomain=S.UniversalSet,
        name='f', **kwargs
    ):
        parameters = Tuple(*[_sympify(a) for a in parameters])

        domain, codomain = _sympify(domain), _sympify(codomain)

        if not isinstance(name, Str):
            name = Str(name)

        obj = super().__new__(cls, parameters, domain, codomain, name)
        return obj

    @property
    def parameters(self):
        return self.args[0]

    @property
    def domain(self):
        return self.args[1]

    @property
    def codomain(self):
        return self.args[2]

    @property
    def name(self):
        return self.args[3]

    @property
    def nargs(self):
        if isinstance(self.domain, ProductSet):
            return len(self.domain.args)
        else:
            return 1

    def eval(self, *args):
        return None

    def __call__(self, *args, evaluate=False, **kwargs):
        if evaluate:
            ret = self.eval(*args)
            if ret is not None:
                return ret
        return AppliedMap(self, *args, evaluate=False)

    def inverse(self, evaluate=False):
        return InverseMap(self, evaluate=evaluate)
    inv = inverse

    def _eval_inverse(self):
        return None

    def doit(self, **hints):
        deep = hints.get('deep', True)
        if deep:
            args = (a.doit(**hints) for a in self.args)
            return self.func(*args, evaluate=True)
        else:
            return self.func(*self.args, evaluate=True)

class InverseMap(Map):
    """
    A class for unevaluated inverse mappings

    Parameters
    ==========

    mapping : Map

    evaluate : bool, optional
        If True, return the evaluated inverse function of *mapping*

    Example
    =======

    >>> from sympy.map import Map, InverseMap
    >>> from sympy.abc import x

    # Minimalist implementation of Exp and Log

    >>> class Exp(Map):
    ...     def _eval_inverse(self):
    ...         return Log(self.parameters)
    >>> exp = Exp(parameters=(2,))

    >>> class Log(Map):
    ...     def _eval_inverse(self):
    ...         return Exp(self.parameters)
    >>> log = Log(parameters=(2,))

    >>> exp.inv() == InverseMap(exp)
    True
    >>> InverseMap(exp, evaluate=True) == exp.inv(evaluate=True) == log
    True

    """
    def __new__(cls, mapping, evaluate=False, **kwargs):
        if evaluate:
            obj = mapping._eval_inverse()
            if obj is None:
                obj = cls(mapping, evaluate=False)
        else:
            obj = super(Expr, cls).__new__(cls, mapping)
        return obj

    @property
    def base(self):
        return self.args[0]

    @property
    def parameters(self):
        return self.base.parameters

    @property
    def domain(self):
        return self.base.codomain

    @property
    def codomain(self):
        return self.base.domain

    def _eval_inverse(self):
        return self.base

class IdentityMap(Map):
    """
    General identity mapping

    Explanation
    ===========

    Identity map is a map that always return its argument.
    For common linear identity function, use LinearIdentityMap.

    Examples
    ========

    >>> from sympy.abc import x, y
    >>> from sympy.map import IdentityMap

    >>> I = IdentityMap()
    >>> I(*(x, y), evaluate=True)
    (x, y)

    """
    def __new__(cls, domain=S.UniversalSet, name='', **kwargs):
        if not isinstance(name, Str):
            name = Str(name)
        return super(Expr, cls).__new__(cls, domain, name)

    @property
    def parameters(self):
        return ()

    @property
    def domain(self):
        return self.args[0]

    @property
    def codomain(self):
        return self.domain

    @property
    def name(self):
        return self.args[1]

    def eval(self, *args):
        return args

    def _eval_inverse(self):
        return self

class AppliedMap(Expr):
    """
    Unevaluated result of Map applied to arguments

    Parameters
    ==========

    map : Map

    args
        arguments applied to *map*

    evaluate : bool, optional
        If True, returns evaluated application of *map* to *args*

    Examples
    ========

    >>> from sympy import symbols
    >>> from sympy.map import Map, AppliedMap
    >>> x, y = symbols('x y')

    >>> class F(Map):
    ...     def eval(self, x):
    ...         a, = self.parameters
    ...         return a*x

    >>> f1 = F(parameters=(x,), name='f')
    >>> isinstance(f1(y), AppliedMap)
    True
    >>> f1(y, evaluate=True)
    x*y

    """
    def __new__(cls, map, *args, evaluate=False, **kwargs):
        args = (_sympify(a) for a in args)
        if evaluate:
            return map(*args, evaluate=True)
        return super().__new__(cls, map, *args)

    @property
    def map(self):
        return self.args[0]

    @property
    def arguments(self):
        return self.args[1:]

    @property
    def parameters(self):
        return self.map.parameters

    @property
    def free_symbols(self):
        result = set()
        result.update(self.parameters.free_symbols)
        for a in self.arguments:
            result.update(a.free_symbols)
        return result

    def doit(self, **hints):
        deep = hints.get('deep', True)
        if deep:
            args = (a.doit(**hints) for a in self.args)
            return self.func(*args, evaluate=True)
        else:
            return self.func(*self.args, evaluate=True)
