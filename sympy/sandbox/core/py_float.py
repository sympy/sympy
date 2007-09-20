from utils import memoizer_immutable_args
from basic import Basic, sympify
from number import Real

class Float(Real, float):

    @memoizer_immutable_args('Float.__new__')
    def __new__(cls, f):
        if isinstance(f, Basic):
            return f.evalf()
        return float.__new__(cls, f)

    def __cmp__(self, other):
        s,o = float(self), float(other)
        if s==o: return 0
        if s<o: return -1
        return 1

    def as_native(self):
        return float(self)

    # float has __int__, __float__
