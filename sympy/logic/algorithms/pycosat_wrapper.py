from __future__ import annotations
from sympy.assumptions.cnf import EncodedCNF

def pycosat_satisfiable(expr, all_models=False):
    import pycosat
    if not isinstance(expr, EncodedCNF):
        exprs = EncodedCNF()
        exprs.add_prop(expr)
        expr = exprs

    if not expr.data:
        if all_models:
            return iter([{True: True}])
        return {True: True}

    if set() in expr.data:
        if all_models:
            return (f for f in [False])
        return False

    if not all_models:
        r = pycosat.solve(expr.data)
        result = (r != "UNSAT")
        if not result:
            return result
        res = {expr.symbols[abs(lit) - 1]: lit > 0 for lit in r}
        return {True: True} if res == {} else res
    else:
        r = pycosat.itersolve(expr.data)

        def _gen(results):
            satisfiable = False
            try:
                while True:
                    sol = next(results)
                    res = {expr.symbols[abs(lit) - 1]: lit > 0 for lit in sol}
                    yield {True: True} if res == {} else res
                    satisfiable = True
            except StopIteration:
                if not satisfiable:
                    yield False

        return _gen(r)
