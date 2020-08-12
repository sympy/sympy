from sympy import S, Basic, sympify, Tuple, Integer
from sympy.core.function import ArgumentIndexError
from .map import Map, AppliedMap, IdentityMap

__all__ = [
    'DiffOp', 'DerivativeFunction',
]

class DiffOp(Map):
    """
    The rudimentary implementation of differential operator.
    After function space, topology, etc are implemented as SymPy object,
    this class will be enhanced.

    """

    # Will be changed after function space is introduced as SymPy object
    domain = codomain = S.UniversalSet

    def __new__(cls, *indices, **kwargs):

        if not indices:
            return IdentityMap(cls.domain)

        # standardizing the indices
        array_likes = (tuple, list, Tuple)
        if isinstance(indices, array_likes) and\
            len(indices) == 1 and isinstance(indices[0], array_likes) and\
            all(isinstance(i, Tuple) for i in indices[0]):
            # this is the case where expr.func(*expr.args) is run
            indices_count = indices[0]
        else:
            indices = list(sympify(indices))
            indices_count = []
            for i,v in enumerate(indices):
                if isinstance(v, Integer):
                    tup = Tuple(v, S.One)
                elif isinstance(v, array_likes):
                    tup = Tuple(*v)
                else:
                    raise TypeError("%s cannot be index of differential operator." % v)

                if i == 0:
                    indices_count.append(tup)
                else:
                    prev, prevcount = indices_count[-1]
                    if prev == tup[0]:
                        indices_count[-1] = Tuple(prev, prevcount+tup[1])
                    else:
                        indices_count.append(tup)
            indices_count = Tuple(*indices_count)

        return super().__new__(cls, indices_count)

    @property
    def indices(self):
        return self.args[0]

    @property
    def free_symbols(self):
        return {}

    def __call__(self, f, **kwargs):
        return DerivativeFunction(self, (f,), **kwargs)

    def apply(self, f, **kwargs):
        nargs = f.nargs
        for index, _ in self.indices:
            if index not in range(1, nargs+1):
                raise ArgumentIndexError(f, 1)
        return super().apply(f, **kwargs)

    def eval(self, f, **kwargs):
        diff_f = f
        for index, count in self.indices:
            for i in range(count):
                diff_f = diff_f.fdiff(index)
                if not isinstance(diff_f, Basic):
                    return None
        return diff_f

    def _eval_subs(self, old, new):
        # Do not allow substitution
        return self

class DerivativeFunction(Map, AppliedMap):

    def _corresponding_oldfunc(self):
        return super()._corresponding_oldfunc() +\
            [Derivative]

    def __new__(cls, map, args, **kwargs):
        return super(Map, DerivativeFunction).__new__(cls, map, args, **kwargs)

    @property
    def operator(self):
        return self.map

    @property
    def function(self):
        return self.arguments[0]

    @property
    def domain(self):
        return self.function.domain

    @property
    def codomain(self):
        return self.function.codomain

    def eval(self, *args, **kwargs):
        fdiff = self.doit(deep=False)
        if fdiff != self:
            return fdiff(*args, evaluate=True)

        # Since we don't know what symbol corresponds to
        # nth index, we can't call Derivative to do the job.
        # (This can be avoided if coordinate system from
        # diffgeom module is integrated with Map.)
        raise NotImplementedError

### imported for Function.__instancecheck__

from sympy.core.function import Derivative
