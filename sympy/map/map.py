from sympy.core import Expr, S, Tuple
from sympy.core.sympify import _sympify
from sympy.core.symbol import Str

class Map(Expr):
    """
    A class for general mathematical mappings

    Parameters
    ==========

    parameters
        Parameters of the map

    domain, codomain : Set

    name : str
        Name of the map

    """
    def __new__(
        cls, parameters, domain=S.UniversalSet, codomain=S.UniversalSet,
        name='', **kwargs
    ):
        parameters = Tuple(*[_sympify(a) for a in parameters])
        domain, codomain = _sympify(domain), _sympify(codomain)
        if not isinstance(name, Str):
            name = Str(name)
        return super().__new__(cls, parameters, domain, codomain, name)

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
