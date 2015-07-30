from sympy.tensor import Idx, IndexedBase, Indexed
from sympy import Integral, Dummy, sympify, Tuple


class _IndexedIntegral(Integral):
    """
    Experimental class to test integration by indexed variables.
    """

    def __new__(cls, function, *limits, **assumptions):
        repl, limits = _IndexedIntegral._indexed_process_limits(limits)
        function = sympify(function)
        function = function.xreplace(repl)
        obj = Integral.__new__(cls, function, *limits, **assumptions)
        obj._indexed_repl = repl
        obj._indexed_reverse_repl = dict((val, key) for key, val in repl.items())
        return obj

    def doit(self):
        res = super(_IndexedIntegral, self).doit()
        return res.xreplace(self._indexed_reverse_repl)

    @staticmethod
    def _indexed_process_limits(limits):
        repl = {}
        newlimits = []
        for i in limits:
            if isinstance(i, (tuple, list, Tuple)):
                v = i[0]
                vrest = i[1:]
            else:
                v = i
                vrest = ()
            if isinstance(v, Indexed):
                if v not in repl:
                    r = Dummy(str(v))
                    repl[v] = r
                newlimits.append((r,)+vrest)
            else:
                newlimits.append(i)
        return repl, newlimits
