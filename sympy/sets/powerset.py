from __future__ import print_function, division

from sympy.core.decorators import _sympifyit
from sympy.core.evaluate import global_evaluate
from sympy.core.logic import fuzzy_bool
from sympy.core.singleton import S
from sympy.core.sympify import _sympify

from .sets import Set, tfn


class PowerSet(Set):
    def __new__(cls, arg, evaluate=global_evaluate[0]):
        arg = _sympify(arg)

        if not isinstance(arg, Set):
            raise ValueError('{} must be a set.'.format(arg))

        if evaluate:
            ret = arg._eval_powerset()

            if ret is not None and not isinstance(arg, PowerSet):
                return ret

        return super(PowerSet, cls).__new__(cls, arg)

    @property
    def arg(self):
        return self.args[0]

    @_sympifyit('other', NotImplemented)
    def _contains(self, other):
        if not isinstance(other, Set):
            raise ValueError('{} must be a set.'.format(other))

        arg = self.arg
        ret = fuzzy_bool(arg.is_superset(other))
        if ret is not None:
            return ret
        return None

    def __len__(self):
        return 2 ** len(self.arg)
