from sympy.printing.smtlib import smtlib_code
from sympy.assumptions.assume import AppliedPredicate
from sympy.assumptions.cnf import EncodedCNF
from sympy.assumptions.ask import Q

from sympy.core import Add, Mul
from sympy.core.relational import Equality, LessThan, GreaterThan, StrictLessThan, StrictGreaterThan
from sympy.core.function import Function, UndefinedFunction, AppliedUndef
from sympy.functions.elementary.complexes import Abs
from sympy.functions.elementary.exponential import Pow
from sympy.functions.elementary.miscellaneous import Min, Max
from sympy.logic.boolalg import And, Or, Xor, Implies
from sympy.logic.boolalg import Not, ITE
from sympy.assumptions.relation.equality import StrictGreaterThanPredicate, StrictLessThanPredicate, GreaterThanPredicate, LessThanPredicate, EqualityPredicate
from sympy.external import import_module
from typing import Callable

# Define supported predicates as a constant set
supported_predicates = {
    Q.gt, Q.lt, Q.ge, Q.le, Q.ne, Q.eq, Q.positive, Q.negative, 
    Q.extended_negative, Q.extended_positive, Q.zero, Q.nonzero, 
    Q.nonnegative, Q.nonpositive, Q.extended_nonzero, 
    Q.extended_nonnegative, Q.extended_nonpositive
}

def is_supported_predicate(pred):
    if not isinstance(pred, AppliedPredicate):
        return False
    
    if pred.function in supported_predicates:
        return True

    return False


def z3_satisfiable(expr, all_models=False, return_model=None):

    if not isinstance(expr, EncodedCNF):
        exprs = EncodedCNF()
        exprs.add_prop(expr)
        expr = exprs

    z3 = import_module("z3")
    if z3 is None:
        raise ImportError("z3 is not installed")
    try:
        s = encoded_cnf_to_z3_solver(expr, z3)
    except Exception:
        return False

    res = str(s.check())
    if res == "unsat":
        return False
    elif res == "sat":
        model = z3_model_to_sympy_model(s.model(), expr)

        has_theory_preds = any(is_supported_predicate(pred) for pred in expr.encoding.keys())

        if return_model is True:
            return model if model else True
        
        elif return_model is False:
            return True
        
        else:
            # - If the CNF has theory predicates => boolean
            # - Else (pure propositional) => model dict (fallback to True if empty)
            if has_theory_preds:
                return True
            else:
                return model if model else True
    else:
        return False

def z3_model_to_sympy_model(z3_model, enc_cnf):
    """Convert Z3 model back to SymPy model."""
    rev_enc = {value: key for key, value in enc_cnf.encoding.items()}
    result = {}
    
    for var in z3_model:
        var_name = var.name()
        # Only process variables that start with 'd' followed by digits
        if var_name.startswith('d') and var_name[1:].isdigit():
            var_id = int(var_name[1:])
            # Only include if this variable ID exists in our reverse encoding
            if var_id in rev_enc:
                result[rev_enc[var_id]] = bool(z3_model[var])
    return result


# Converts a single CNF clause to SMT-LIB assertion format.
def clause_to_assertion(clause):
    clause_strings = [f"d{abs(lit)}" if lit > 0 else f"(not d{abs(lit)})" for lit in clause]
    return "(assert (or " + " ".join(clause_strings) + "))"

def find_uninterpreted_functions(expr):
    uf_arity = {}

    # Handle bare function heads (symbols that are functions but not called)
    for head in expr.atoms(UndefinedFunction):
        uf_arity.setdefault(str(head), None)

    # Handle function calls like f(a), g(x,y)
    if AppliedUndef is not None:
        calls = expr.atoms(AppliedUndef)
    else:
        calls = [c for c in expr.atoms(Function) if isinstance(c.func, UndefinedFunction)]

    for call in calls:
        name = str(call.func)
        arity = len(call.args)
        prev = uf_arity.get(name)
        uf_arity[name] = arity if prev is None else max(prev or 0, arity)

    return uf_arity


def collect_all_uninterpreted_functions(enc_cnf):
    uf_table = {}
    
    for pred, enc in enc_cnf.encoding.items():
        if not is_supported_predicate(pred):
            continue
            
        pred_ufs = find_uninterpreted_functions(pred)
        # Merge with existing, taking maximum arity for each function
        for name, arity in pred_ufs.items():
            if name in uf_table:
                if uf_table[name] is None or arity is None:
                    uf_table[name] = max(uf_table[name] or 0, arity or 0)
                else:
                    uf_table[name] = max(uf_table[name], arity)
            else:
                uf_table[name] = arity
                
    return uf_table


def collect_all_free_symbols(enc_cnf):
    """Collect all free symbols from all predicates in the CNF."""
    free_symbols = set()
    
    for pred, enc in enc_cnf.encoding.items():
        if not is_supported_predicate(pred):
            continue
            
        free_symbols |= pred.free_symbols
        
    return free_symbols


def encoded_cnf_to_z3_solver(enc_cnf, z3):
    s = z3.Solver()
    
    uf_table = collect_all_uninterpreted_functions(enc_cnf)
    free_symbols = collect_all_free_symbols(enc_cnf)
    
    declarations = []
    assertions = []
    
    declarations.extend(f"(declare-const d{var} Bool)" for var in enc_cnf.variables)
    assertions.extend(clause_to_assertion(clause) for clause in enc_cnf.data)
    
    for name, arity in sorted(uf_table.items()):
        if arity is None:
            continue
        elif arity == 0:
            declarations.append(f"(declare-const {name} Real)")
        else:
            arg_types = " ".join(["Real"] * arity)
            declarations.append(f"(declare-fun {name} ({arg_types}) Real)")
    
    for sym in sorted(free_symbols, key=str):
        declarations.append(f"(declare-const {sym} Real)")
    
    for pred, enc in enc_cnf.encoding.items():
        if not is_supported_predicate(pred):
            continue
        
        symbol_table = {}
        
        for sym in pred.free_symbols:
            symbol_table[sym] = float
        
        pred_ufs = find_uninterpreted_functions(pred)
        for name, arity in pred_ufs.items():
            if arity is None or arity == 0:
                symbol_table[Function(name)] = float
            elif arity == 1:
                symbol_table[Function(name)] = Callable[[float], float]
            elif arity == 2:
                symbol_table[Function(name)] = Callable[[float, float], float]
            else:
                symbol_table[Function(name)] = Callable[..., float]
        
        pred_str = smtlib_code(
            pred, 
            auto_declare=False, 
            auto_assert=False, 
            known_functions=known_functions,
            symbol_table=symbol_table
        )
        
        clause = f"(implies d{enc} {pred_str})"
        assertion = f"(assert {clause})"
        assertions.append(assertion)
    
    declarations = "\n".join(declarations)
    assertions= "\n".join(assertions)
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
