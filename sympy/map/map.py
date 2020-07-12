from sympy.core import Expr, S, Tuple
from sympy.core.compatibility import is_sequence, ordered
from sympy.core.function import arity
from sympy.core.sympify import _sympify
from sympy.core.symbol import Str
from sympy.sets.sets import FiniteSet

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
    >>> x, y, z = symbols('x y z')

    >>> f = Map(parameters=(x, ), name='f')
    >>> f.name
    'f'
    >>> f.parameters
    (x,)
    >>> f.nargs
    Naturals0

    >>> class G(Map):
    ...     def eval(self, a, b):
    ...         x,y = self.parameters
    ...         return a*x + b*y
    >>> g = G(parameters=(x,y), name='g')

    >>> g.name
    'g'
    >>> g.parameters
    (x, y)
    >>> g.nargs
    FiniteSet(2)

    >>> g.subs(x,2).parameters
    (2, y)
    >>> g.subs(y, 4).parameters
    (x, 4)

    """
    def __new__(
        cls, parameters=(), domain=S.UniversalSet, codomain=S.UniversalSet,
        name='f', **kwargs
    ):
        parameters = Tuple(*[_sympify(a) for a in parameters])

        domain, codomain = _sympify(domain), _sympify(codomain)

        if not isinstance(name, Str):
            name = Str(name)

        nargs = cls._find_nargs()

        obj = super().__new__(cls, parameters, domain, codomain, name)
        obj._nargs = nargs
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

    @classmethod
    def _find_nargs(cls):
        nargs = arity(cls.eval)
        if is_sequence(nargs): # multiple arity
            nargs = [i-1 for i in nargs]
            nargs = tuple(ordered(nargs))
        elif nargs is not None:
            nargs = (nargs-1,)
        return FiniteSet(*nargs) if nargs else S.Naturals0

    @property
    def nargs(self):
        return self._nargs

    def eval(self, *args):
        return None

    def __call__(self, *args, evaluate=False, **kwargs):
        if evaluate:
            ret = self.eval(*args)
            if ret is not None:
                return ret
        return AppliedMap(self, *args, evaluate=False)

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
    >>> from sympy.map import Map
    >>> x, y, z = symbols('x y z')

    >>> f = Map(parameters=(x, y))
    >>> expr = f.subs(x, 4)(z)
    >>> expr.arguments
    (z,)
    >>> expr.parameters
    (4, y)
    >>> expr.doit() == expr
    True

    >>> class G1(Map):
    ...     def eval(self, x, y):
    ...         a, b = self.parameters
    ...         return a*x + b*y
    >>> g1 = G1(parameters=(1, 2))
    >>> g1(x, y) != g1(x, y, evaluate=True) == g1(x, y).doit() == x+2*y
    True

    >>> class G2(Map):
    ...     def eval(self, x, y):
    ...         return None
    >>> g2 = G2(parameters=(1, 2))
    >>> g2(x, y) == g2(x, y, evaluate=True) == g2(x, y).doit()
    True

    """
    def __new__(cls, map, *args, evaluate=False):
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
            arguments = (a.doit(**hints) for a in self.arguments)
            return self.func(self.map, *arguments, evaluate=True)
        else:
            return self.func(self.map, *self.arguments, evaluate=True)
