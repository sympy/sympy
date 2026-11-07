from __future__ import annotations
from sympy.assumptions.cnf import EncodedCNF

def minisat22_satisfiable(expr, all_models=False, minimal=False):

    if not isinstance(expr, EncodedCNF):
        exprs = EncodedCNF()
        exprs.add_prop(expr)
        expr = exprs

    from pysat.solvers import Minisat22

    # Return UNSAT when False (encoded as 0) or empty clause is present in the CNF
    if {0} in expr.data or set() in expr.data:
        if all_models:
            return (f for f in [False])
        return False

    r = Minisat22(expr.data)

    if minimal:
        r.set_phases([-(i+1) for i in range(r.nof_vars())])

    if not all_models:
        if not r.solve():
            return False
        res = {expr.symbols[abs(lit) - 1]: lit > 0 for lit in r.get_model()}
        return {True: True} if res == {} else res

    else:
        # Make solutions SymPy compatible by creating a generator
        def _gen(results):
            satisfiable = False
            while results.solve():
                sol = results.get_model()
                res = {expr.symbols[abs(lit) - 1]: lit > 0 for lit in sol}
                yield {True: True} if res == {} else res
                if minimal:
                    results.add_clause([-i for i in sol if i>0])
                else:
                    results.add_clause([-i for i in sol])
                satisfiable = True
            if not satisfiable:
                yield False

        return _gen(r)
