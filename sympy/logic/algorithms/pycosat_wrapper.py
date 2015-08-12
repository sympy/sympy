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
        r = pycosat.solve(map(list, clauses_int_repr))
        result = (r != "UNSAT")
        if not result:
            return result
        if expr:
            return dict((symbols[abs(lit) - 1], lit > 0) for lit in r)
        return bool(expr)
    else:
        r = pycosat.itersolve(map(list, clauses_int_repr))
        result = (r != "UNSAT")
        if not result:
            return result
        sols = []
        # Make solutions sympy compatible
        for sol in r:
            sols.append(dict((symbols[abs(lit) - 1], lit > 0) for lit in sol))
        # Return generator object
        if expr:
            return (sol for sol in sols)
        return (i for i in [bool(expr)])
