
from basic import Basic
from sympify import _sympify

class Interval(Basic):

    def __new__(cls, start, end, **assumptions):
        start = _sympify(start)
        end = _sympify(end)
        return Basic.__new__(cls, start, end, **assumptions)

    @property
    def start(self):
        return self._args[0]

    @property
    def end(self):
        return self._args[1]

# /cyclic/
import sympify as _
_.Interval  = Interval
del _
