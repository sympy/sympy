from sympy.printing.smtlib import smtlib_code
from sympy import symbols
from sympy.core.function import Function, UndefinedFunction, AppliedUndef
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

n = 4
uninterpreted_functions = symbols(f"f1:{n+1}", cls=Function)   # (f1, f2, f3, f4)
uninterpreted_variables = symbols(f"u1:{n+1}")   # (u1, u2, u3, u4)

x = symbols('x')

# The printer spits out the smt lib code, but WE disabled the auto_declare fucntionality

def find_uninterpreted_functions(expr):
    # takes all atoms, spits out the symbol table for everything
    # f: []
    uf_arity = {}
    # bare heads, e.g. a symbol 'f' that shows up without being called
    for head in expr.atoms(UndefinedFunction):
        uf_arity[head] = Callable[[float], float]

    return uf_arity


def encoded_cnf_to_z3_solver(enc_cnf, z3):
    def dummify_bool(pred):
        return False
        assert isinstance(pred, AppliedPredicate)

        if pred.function in [Q.positive, Q.negative, Q.zero]:
            return pred
        else:
            return False

    s = z3.Solver()

    uf_table = {}
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


        uf_table = find_uninterpreted_functions(pred)
        print(f"uf_table = {uf_table}")
        pred_str = smtlib_code(pred, auto_declare=False, known_types=known_types, auto_assert=False, known_functions=known_functions, symbol_table=uf_table)

        free_symbols |= pred.free_symbols
        pred = pred_str
        clause = f"(implies d{enc} {pred})"
        assertion = "(assert " + clause + ")"
        assertions.append(assertion)

    # assertions.append(f"(assert (implies d{enc} {pred_str}))")


    for sym in free_symbols:
        declarations.append(f"(declare-const {sym} Real)")

    # declare uninterpreted functions for real numbers
    for sym in uf_table.keys():
        declarations.append(f"declare-fun {sym} (Real) Real")
    

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