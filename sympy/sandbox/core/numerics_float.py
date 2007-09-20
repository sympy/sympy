from utils import memoizer_immutable_args
from basic import Basic, sympify
from number import Real

from numerics import Float as num_Float

class Float(Real, num_Float):

    @memoizer_immutable_args('Float.__new__')
    def __new__(cls, f):
        if isinstance(f, Basic):
            return f.evalf()
        return num_Float.__new__(cls, f)

    def as_native(self):
        return num_Float(self)

