from sympy.printing.smtlib import smtlib_code
from sympy import symbols
from sympy.core.function import Function
from sympy.assumptions.assume import AppliedPredicate
from sympy.assumptions.cnf import EncodedCNF
from sympy.assumptions.ask import Q

from sympy.core import Add, Mul, Function, Number
from sympy.core.relational import Equality, LessThan, GreaterThan, StrictLessThan, StrictGreaterThan
from sympy.functions.elementary.complexes import Abs
from sympy.functions.elementary.exponential import Pow
from sympy.functions.elementary.miscellaneous import Min, Max
from sympy.logic.boolalg import And, Or, Xor, Implies
from sympy.logic.boolalg import Not, ITE
from sympy.assumptions.relation.equality import StrictGreaterThanPredicate, StrictLessThanPredicate, GreaterThanPredicate, LessThanPredicate, EqualityPredicate
from sympy.external import import_module
from typing import Callable

def z3_satisfiable(expr, all_models=False):
    if not isinstance(expr, EncodedCNF):
        exprs = EncodedCNF()
        exprs.add_prop(expr)
        expr = exprs

    z3 = import_module("z3")
    if z3 is None:
        raise ImportError("z3 is not installed")

    s = encoded_cnf_to_z3_solver(expr, z3)

    res = str(s.check())
    if res == "unsat":
        return False
    elif res == "sat":
        return z3_model_to_sympy_model(s.model(), expr)
    else:
        return None


def z3_model_to_sympy_model(z3_model, enc_cnf):
    rev_enc = {value : key for key, value in enc_cnf.encoding.items()}
    return {rev_enc[int(var.name()[1:])] : bool(z3_model[var]) for var in z3_model}


def clause_to_assertion(clause):
    clause_strings = [f"d{abs(lit)}" if lit > 0 else f"(not d{abs(lit)})" for lit in clause]
    return "(assert (or " + " ".join(clause_strings) + "))"


def _infer_dom_sort(enc_cnf) -> str:
    """Return 'Real' if numeric relations/terms appear; else 'U'."""
    numeric_rels = {
        Q.gt, Q.ge, Q.lt, Q.le,
        Q.positive, Q.negative, Q.nonnegative, Q.nonpositive,
        Q.extended_positive, Q.extended_negative,
        Q.extended_nonnegative, Q.extended_nonpositive,
    }
    for pred in enc_cnf.encoding.keys():
        if isinstance(pred, AppliedPredicate) and getattr(pred, "function", None) in numeric_rels:
            return "Real"
        # Any numeric structure/literals inside the predicate?
        if pred.atoms(Number) or pred.atoms(Add) or pred.atoms(Mul) or pred.atoms(Pow):
            return "Real"
    return "U"

from sympy.printing.smtlib import smtlib_code
from sympy.assumptions.assume import AppliedPredicate
from sympy.assumptions.cnf import EncodedCNF
from sympy.assumptions.ask import Q

from sympy.core import Add, Mul
from sympy.core.relational import Equality, LessThan, GreaterThan, StrictLessThan, StrictGreaterThan
from sympy.functions.elementary.complexes import Abs
from sympy.functions.elementary.exponential import Pow
from sympy.functions.elementary.miscellaneous import Min, Max
from sympy.logic.boolalg import And, Or, Xor, Implies
from sympy.logic.boolalg import Not, ITE
from sympy.assumptions.relation.equality import StrictGreaterThanPredicate, StrictLessThanPredicate, GreaterThanPredicate, LessThanPredicate, EqualityPredicate
from sympy.external import import_module

def z3_satisfiable(expr, all_models=False):
    if not isinstance(expr, EncodedCNF):
        exprs = EncodedCNF()
        exprs.add_prop(expr)
        expr = exprs

    z3 = import_module("z3")
    if z3 is None:
        raise ImportError("z3 is not installed")

    s = encoded_cnf_to_z3_solver(expr, z3)

    res = str(s.check())
    if res == "unsat":
        return False
    elif res == "sat":
        return z3_model_to_sympy_model(s.model(), expr)
    else:
        return None


def z3_model_to_sympy_model(z3_model, enc_cnf):
    rev_enc = {value : key for key, value in enc_cnf.encoding.items()}
    return {rev_enc[int(var.name()[1:])] : bool(z3_model[var]) for var in z3_model}


def clause_to_assertion(clause):
    clause_strings = [f"d{abs(lit)}" if lit > 0 else f"(not d{abs(lit)})" for lit in clause]
    return "(assert (or " + " ".join(clause_strings) + "))"


def encoded_cnf_to_z3_solver(enc_cnf, z3):
    def dummify_bool(pred):
        return False
        assert isinstance(pred, AppliedPredicate)

        if pred.function in [Q.positive, Q.negative, Q.zero]:
            return pred
        else:
            return False

    n = 4
    uninterpreted_functions = symbols(f"f1:{n+1}", cls=Function)   # (f1, f2, f3, f4)
    x = symbols('x')

    class U: 
        pass

    symbol_table = {x: U}
    for f in uninterpreted_functions:
        symbol_table[f] = Callable[[U], U]

    known_types = {U: "U"}

    s = z3.Solver()

    declarations = []
    declarations.append(f"(declare-sort U 0)")
    for fhead in uninterpreted_functions:
        declarations.append(f"(declare-fun {fhead} (U) U)")

    declarations += [f"(declare-const d{var} Bool)" for var in enc_cnf.variables]
    assertions = [clause_to_assertion(clause) for clause in enc_cnf.data]

    # track constants by sort
    uf_arg_consts = set()   # symbols used as arguments to any UF -> declare as U
    real_consts   = set()   # other free symbols -> declare as Real
    free_symbols = set()
    for pred, enc in enc_cnf.encoding.items():
        if not isinstance(pred, AppliedPredicate):
            continue
        if pred.function not in (Q.gt, Q.lt, Q.ge, Q.le, Q.ne, Q.eq, Q.positive, Q.negative, Q.extended_negative, Q.extended_positive, Q.zero, Q.nonzero, Q.nonnegative, Q.nonpositive, Q.extended_nonzero, Q.extended_nonnegative, Q.extended_nonpositive):
            continue

        args_used_by_ufs = set()
        for call in pred.atoms(Function):
            if call.func in uninterpreted_functions:
                for arg in call.args:
                    if getattr(arg, "is_Symbol", False):
                        args_used_by_ufs.add(arg)
        # update global sets
        uf_arg_consts |= args_used_by_ufs
        real_consts   |= (pred.free_symbols - args_used_by_ufs)

        # local typing for this predicate: UF args -> U, others -> Real
        symbol_table = dict(symbol_table)
        for sym in args_used_by_ufs:
            symbol_table[sym] = U
        for sym in (pred.free_symbols - args_used_by_ufs):
            symbol_table[sym] = "Real"

        pred_str = smtlib_code(pred, auto_declare=False, known_types=known_types, auto_assert=False, known_functions=known_functions, symbol_table=symbol_table)

        free_symbols |= pred.free_symbols
        pred = pred_str
        clause = f"(implies d{enc} {pred})"
        assertion = "(assert " + clause + ")"
        assertions.append(assertion)

    # assertions.append(f"(assert (implies d{enc} {pred_str}))")


    for sym in free_symbols:
        declarations.append(f"(declare-const {sym} Real)")
    

    print(f"declarations: \n{declarations}")
    print(f"assertions: \n{assertions}")

    declarations = "\n".join(declarations)
    assertions = "\n".join(assertions)
    s.from_string(declarations)
    s.from_string(assertions)

    return s


known_functions = {
            Add: '+',
            Mul: '*',

            Equality: '=',
            LessThan: '<=',
            GreaterThan: '>=',
            StrictLessThan: '<',
            StrictGreaterThan: '>',

            EqualityPredicate(): '=',
            LessThanPredicate(): '<=',
            GreaterThanPredicate(): '>=',
            StrictLessThanPredicate(): '<',
            StrictGreaterThanPredicate(): '>',

            Abs: 'abs',
            Min: 'min',
            Max: 'max',
            Pow: '^',

            And: 'and',
            Or: 'or',
            Xor: 'xor',
            Not: 'not',
            ITE: 'ite',
            Implies: '=>',
        }