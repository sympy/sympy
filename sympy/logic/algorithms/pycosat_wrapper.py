from __future__ import print_function, division

from collections import defaultdict
from heapq import heappush, heappop

from sympy import default_sort_key
from sympy.logic.boolalg import conjuncts, to_cnf, to_int_repr, _find_predicates, Boolean


def pycosat_satisfiable(expr, all_models):
    import pycosat
    symbols = sorted(_find_predicates(expr), key=default_sort_key)
    symbols_int_repr = range(1, len(symbols) + 1)
    clauses = conjuncts(to_cnf(expr))
    clauses_int_repr = to_int_repr(clauses, symbols)

    if not all_models:
        r = pycosat.solve(clauses_int_repr)
        result = (r != "UNSAT")
        if not result:
            return result
        if expr:
            return dict((symbols[abs(lit) - 1], lit > 0) for lit in r)
        return bool(expr)
    else:
        r = pycosat.itersolve(clauses_int_repr)
        result = (r != "UNSAT")
        if not result:
            return result

        # Make solutions sympy compatible by creating a generator
        def _gen(results):
          return (dict((symbols[abs(lit) - 1], lit > 0) for lit in sol) for sol in results)

        if expr:
            return _gen(r)
        return (i for i in [bool(expr)])
