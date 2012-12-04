from .basic import Basic
from .singleton import S
from .sympify import _sympify
from .evaluate import global_evaluate
from .expr import Expr


class NewFunction(Basic):
    def __new__(cls,
        *args, domain=S.UniversalSet, codomain=S.UniversalSet):
        domain = _sympify(domain)
        codomain = _sympify(codomain)
        args = [_sympify(arg) for arg in args]
        return Basic.__new__(cls, *args, domain, codomain)

    @property
    def domain(self):
        return self.args[-2]

    @property
    def codomain(self):
        return self.args[-1]

    def eval(self, *args):
        return None

    def __call__(self, *args, evaluate=global_evaluate[0]):
        args = [_sympify(arg) for arg in args]
        if evaluate:
            ret = self.eval(*args)
            if ret is not None:
                return ret
        return AppliedFunction(self, *args, evaluate=False)


class AppliedFunction(Expr):
    def __new__(cls, func, *args, evaluate=global_evaluate[0]):
        if evaluate:
            return func.__call__(*args)
        return Expr.__new__(cls, func, *args)
