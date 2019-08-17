from __future__ import print_function, division

from collections import defaultdict
from heapq import heappush, heappop

from sympy import default_sort_key
from sympy.assumptions.cnf import EncodedCNF
from sympy.logic.boolalg import conjuncts, to_cnf, to_int_repr, _find_predicates, Boolean


def pycosat_satisfiable(expr, all_models=False):
    import pycosat
    if not isinstance(expr, EncodedCNF):
        exprs = EncodedCNF()
        exprs.add_prop(expr)
        expr = exprs

    # Return UNSAT when False (encoded as 0) is present in the CNF
    if {0} in expr.data:
        if all_models:
            return (f for f in [False])
        return False

    if not all_models:
        r = pycosat.solve(expr.data)
        result = (r != "UNSAT")
        if not result:
            return result
        if expr:
            return dict((expr.symbols[abs(lit) - 1], lit > 0) for lit in r)
        return bool(expr)
    else:
        r = pycosat.itersolve(expr.data)
        result = (r != "UNSAT")
        if not result:
            return result

        # Make solutions sympy compatible by creating a generator
        def _gen(results):
          return (dict((expr.symbols[abs(lit) - 1], lit > 0) for lit in sol) for sol in results)

        if expr:
            return _gen(r)
        return (i for i in [bool(expr)])
