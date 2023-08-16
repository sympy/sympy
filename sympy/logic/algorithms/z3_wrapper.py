from sympy.printing.smtlib import smtlib_code
from sympy.assumptions.assume import AppliedPredicate
from sympy.assumptions.cnf import EncodedCNF
def z3_satisfiable(expr, all_models=False):
    if not isinstance(expr, EncodedCNF):
        exprs = EncodedCNF()
        exprs.add_prop(expr)
        expr = exprs

    import z3

    s = encoded_cnf_to_z3_solver(expr, z3)

    res = str(s.check())
    if res == "unsat":
        return False
    elif res == "sat":
        return True
    else:
        return None




def encoded_cnf_to_z3_solver(enc_cnf, z3):
    def dummify_bool(pred):
        return False
        assert isinstance(pred, AppliedPredicate)

        if pred.function in [Q.positive, Q.negative, Q.zero]:
            return pred
        else:
            return False

    s = z3.Solver()

    declarations = []
    for var in enc_cnf.variables:
        declarations.append(f"(declare-const d{var} Bool)")

    declarations = "\n".join(declarations)

    assertions = []
    for clause in enc_cnf.data:
        clause_strings = []
        for lit in clause:
            if lit > 0:
                clause_strings.append(f"d{abs(lit)}")
            elif lit < 0:
                clause_strings.append(f"(not d{abs(lit)})")
            else:
                assert False

        clause = " ".join(clause_strings)
        clause = "(or " + clause + ")"
        assertion = "(assert " + clause + ")"
        assertions.append(assertion)

    assertions = "\n".join(assertions)
    #print(declarations + "\n" + assertions)
    s.from_string(declarations + assertions)

    return s
