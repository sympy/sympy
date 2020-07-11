from sympy.core import Expr, S
from sympy.core.sympify import _sympify
from sympy.core.symbol import Str

class Map(Expr):
    """
    A class for general mathematical mappings
    """
    def __new__(
        cls, *args, domain=S.UniversialSet, codomain=S.UniversialSet,
        name='', **kwargs
    ):
        domain, codomain = _sympify(domain), _sympify(codomain)
        args = (_sympify(a) for a in args)
        name = Str(name)
        return super().__new__(cls, *args, domain, codomain, name)

    @property
    def domain(self):
        return self.args[-3]

    @property
    def codomain(self):
        return self.args[-2]

    @property
    def name(self):
        return self.args[-1]

    def eval(self, *args):
        return None

    def __call__(self, *args, evaluate=False, **kwargs):
        if evaluate:
            ret = self.eval(*args)
            if ret is not None:
                return ret
        return AppliedMap(self, *args, evalaute=False)

class AppliedMap(Expr):
    def __new__(cls, func, *args, evalaute=False):
        if evaluate:
            return func(*args, evalaute=True)
        return super().__new__(cls, func, *args)
