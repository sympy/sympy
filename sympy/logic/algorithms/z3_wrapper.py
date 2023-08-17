from sympy.printing.smtlib import smtlib_code
from sympy.assumptions.assume import AppliedPredicate
from sympy.assumptions.cnf import EncodedCNF
from sympy.assumptions.ask import Q

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

    symbols = set()
    for pred, enc in enc_cnf.encoding.items():
        if not isinstance(pred, AppliedPredicate):
            continue
        if pred.function not in (Q.gt, Q.lt, Q.ge, Q.le, Q.ne, Q.eq, Q.positive, Q.negative, Q.extended_negative, Q.extended_positive, Q.zero, Q.nonzero, Q.nonnegative, Q.nonpositive, Q.extended_nonzero, Q.extended_nonnegative, Q.extended_nonpositive):
            continue


        try:
            pred_str = smtlib_code(pred, auto_declare=False, auto_assert=False)
        except KeyError:
            # Sometimes pred can contain a matrix or something else the smtlib printer
            # doesn't know how to print.
            continue

        symbols |= pred.free_symbols
        pred = pred_str
        clause = f"(implies d{enc} {pred})"
        assertion = "(assert " + clause + ")"
        assertions.append(assertion)

    for sym in symbols:
        declarations.append(f"(declare-const {sym} Real)")


    declarations = "\n".join(declarations)
    assertions = "\n".join(assertions)
    #print(declarations + "\n" + assertions)
    s.from_string(declarations + assertions)

    return s
