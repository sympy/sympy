from __future__ import annotations

from sympy import (
    Symbol, Function, Integer, Float, Eq, Ne, Not, And, Or,
    Implies, Xor, Piecewise, Add, Mul, Mod, StrictLessThan,
    StrictGreaterThan, LessThan, GreaterThan, Expr, S, floor
)
from sympy.assumptions import Q
from sympy.functions.elementary.complexes import Abs
from sympy.external import import_module
from sympy.parsing.smtlib.lark.smtlib_parser import (
    UnknownSMTLibCommandError, UnknownSMTLibOperatorError
)

lark = import_module("lark")

if lark:
    from lark import Transformer, Token, Tree  # type: ignore
else:
    class Transformer:  # type: ignore
        def transform(self, *args):
            pass

    class Token:  # type: ignore
        pass

    class Tree:  # type: ignore
        pass

# Operators from SMT-LIB theories (bitvectors, strings, arrays) that have no
# SymPy representation. These raise NotImplementedError rather than silently
# producing an uninterpreted Function that looks like a regular expression
# but has none of the theory's semantics.
_UNSUPPORTED_THEORY_OPS = frozenset({
    'concat', 'bvnot', 'bvand', 'bvor', 'bvneg', 'bvadd', 'bvmul',
    'bvudiv', 'bvurem', 'bvshl', 'bvlshr', 'bvsub', 'bvult', 'bvxor',
    'bvnand', 'bvnor', 'bvxnor', 'bvcomp', 'bvsdiv', 'bvsrem',
    'bvsmod', 'bvashr', 'bvule', 'bvugt', 'bvuge', 'bvslt', 'bvsle',
    'bvsgt', 'bvsge', 'bv2nat', 'bv2int', 'ubv_to_int',
    'ext_rotate_left', 'ext_rotate_right', 'str.len', 'str.++',
    'str.at', 'str.contains', 'str.indexof', 'str.replace',
    'str.substr', 'str.prefixof', 'str.suffixof', 'str.to.int',
    'int.to.str', 'select', 'store', 'const'
})


class SMTLibTransformer(Transformer):
    """
    Returns a tuple of (symbols, assertions) by evaluating the SMT-LIB AST.
    """

    def __init__(self):
        super().__init__()
        self._symbols = {}
        self._functions = {}
        self._macros = {}
        self._assertions = []
        self._logic = None
        self._sorts = {}
        self._sort_aliases = {}
        self._options = {}
        self._info = {}

        self._assertion_stack = []
        self._symbol_stack = []
        self._macro_stack = []
        self._function_stack = []
        self._sort_stack = []
        self._sort_alias_stack = []

    def get_result(self):
        return self._symbols, self._assertions

    def _resolve_sort(self, sort):
        # Resolve 0-ary define-sort aliases like (define-sort MyInt () Int)
        seen = set()
        while isinstance(sort, str) and sort in self._sort_aliases and sort not in seen:
            seen.add(sort)
            sort = self._sort_aliases[sort]
        return sort

    # Tokens
    def SYMBOL(self, token):
        return str(token)

    def QUOTED_SYMBOL(self, token):
        # |a b| denotes the symbol "a b"; the pipes are quoting syntax, not
        # part of the name
        return str(token)[1:-1]

    def KEYWORD(self, token):
        return str(token)

    def NUMERAL(self, token):
        return Integer(token)

    def DECIMAL(self, token):
        return Float(token)

    def HEXADECIMAL(self, token):
        # TODO: Flagged for future BitVector (BV) work.
        # Currently, bit-width is lost when casting `#x0A` directly to `Integer(10)`.
        return Integer(int(str(token).replace('#x', '0x'), 0))

    def BINARY(self, token):
        # TODO: Flagged for future BitVector (BV) work.
        # Currently, bit-width is lost when casting `#b1010` directly to `Integer(10)`.
        return Integer(int(str(token).replace('#b', '0b'), 0))

    def STRING(self, token):
        # Strip the quotes and undo the "" escape. The result is returned as
        # a Symbol because SymPy doesn't have a native String type.
        return Symbol(str(token)[1:-1].replace('""', '"'))

    # Non-terminals
    def symbol(self, args):
        return args[0]

    def keyword(self, args):
        return args[0]

    def numeral(self, args):
        return args[0]

    def decimal(self, args):
        return args[0]

    def hexadecimal(self, args):
        return args[0]

    def binary(self, args):
        return args[0]

    def string(self, args):
        return args[0]

    def spec_constant(self, args):
        return args[0]

    def sort_simple(self, args):
        return args[0]

    def sort_param(self, args):
        return tuple(args)

    def sort_indexed(self, args):
        # args = [symbol, numeral, ...]
        # We model (_ BitVec 32) as a tuple like parameterised sorts
        return ("_", args[0]) + tuple(args[1:])

    def sorted_var(self, args):
        return (args[0], args[1])

    def var_binding(self, args):
        return (args[0], args[1])

    def list(self, args):
        return args

    # Term handlers
    def term_spec_constant(self, args):
        return args[0]

    def qual_id_simple(self, args):
        return args[0]

    def _indexed_identifier(self, name):
        # (_ bvN width) is a bitvector literal with value N. The width is
        # currently discarded, like for #x and #b literals.
        if name.startswith('bv') and name[2:].isdigit():
            return Integer(name[2:])
        raise NotImplementedError(
            f"Indexed identifier (_ {name} ...) is not supported")

    def qual_id_indexed(self, args):
        return self._indexed_identifier(str(args[0]))

    def qual_id_as(self, args):
        return args[0]

    def qual_id_as_indexed(self, args):
        return self._indexed_identifier(str(args[0]))

    def term_qual_identifier(self, args):
        atom = args[0]
        # Indexed identifiers like (_ bv5 32) have already been converted to
        # SymPy objects
        if isinstance(atom, Expr):
            return atom

        if atom in self._symbols:
            return self._symbols[atom]

        if atom == 'true':
            return S.true
        if atom == 'false':
            return S.false

        return Symbol(str(atom))

    def term_apply(self, args):
        # args[0] is qual_identifier
        op_node = args[0]
        params = args[1:]

        if isinstance(op_node, Expr):
            raise UnknownSMTLibOperatorError(
                f"Cannot apply arguments to expression {op_node}")

        op = str(op_node)

        if op in self._macros:
            macro_args, macro_body = self._macros[op]
            return macro_body.xreplace(dict(zip(macro_args, params)))

        if op in self._functions:
            return self._functions[op](*params)

        # SMT-LIB symbols are case-sensitive, so operators are matched
        # exactly
        if op == 'not':
            return Not(params[0])
        elif op == 'and':
            return And(*params)
        elif op == 'or':
            return Or(*params)
        elif op == '=>':
            result = params[-1]
            for p in reversed(params[:-1]):
                result = Implies(p, result)
            return result
        elif op == 'xor':
            return Xor(*params)
        elif op == '=':
            if len(params) == 2:
                return Eq(params[0], params[1])
            return And(*[Eq(params[i], params[i+1]) for i in range(len(params)-1)])
        elif op == 'distinct':
            if len(params) == 2:
                return Ne(params[0], params[1])
            exprs = []
            for i in range(len(params)):
                for j in range(i + 1, len(params)):
                    exprs.append(Ne(params[i], params[j]))
            return And(*exprs)
        elif op in ('<', '<=', '>', '>='):
            if len(params) == 2:
                if op == '<': return StrictLessThan(params[0], params[1])
                if op == '<=': return LessThan(params[0], params[1])
                if op == '>': return StrictGreaterThan(params[0], params[1])
                if op == '>=': return GreaterThan(params[0], params[1])
            exprs = []
            for i in range(len(params) - 1):
                if op == '<': exprs.append(StrictLessThan(params[i], params[i+1]))
                elif op == '<=': exprs.append(LessThan(params[i], params[i+1]))
                elif op == '>': exprs.append(StrictGreaterThan(params[i], params[i+1]))
                elif op == '>=': exprs.append(GreaterThan(params[i], params[i+1]))
            return And(*exprs)
        elif op == '+':
            return Add(*params)
        elif op == '-':
            if len(params) == 1:
                return -params[0]
            return params[0] - Add(*params[1:])
        elif op == '*':
            return Mul(*params)
        elif op == '/':
            if len(params) == 2:
                return params[0] / params[1]
            return params[0] / Mul(*params[1:])
        elif op == 'div':
            # SMT-LIB integer division is Euclidean: the remainder is always
            # nonnegative, which matches floor division only for positive
            # divisors
            result = params[0]
            for p in params[1:]:
                result = (result - Mod(result, Abs(p))) / p
            return result
        elif op == 'mod':
            return Mod(params[0], Abs(params[1]))
        elif op == 'ite':
            return Piecewise((params[1], params[0]), (params[2], True))
        elif op == 'abs':
            return Abs(params[0])
        elif op == 'to_real':
            # SymPy does not distinguish an integer from the corresponding
            # real number
            return params[0]
        elif op == 'to_int':
            return floor(params[0])
        elif op == 'is_int':
            return Eq(params[0], floor(params[0]))

        if op in _UNSUPPORTED_THEORY_OPS:
            raise NotImplementedError(
                f"The SMT-LIB operator '{op}' has no SymPy equivalent and is "
                "not supported")

        raise UnknownSMTLibOperatorError(f"Unknown operator or function: {op}")

    def term_let(self, args):
        bindings = args[:-1]
        term = args[-1]
        subs_dict = {Symbol(var): val for var, val in bindings}
        return term.xreplace(subs_dict)

    def term_forall(self, args):
        raise NotImplementedError(
            "Quantifiers (forall) are not supported since SymPy has no "
            "representation for them")

    def term_exists(self, args):
        raise NotImplementedError(
            "Quantifiers (exists) are not supported since SymPy has no "
            "representation for them")

    def term_match(self, args):
        raise NotImplementedError("Match expressions are not fully supported yet.")

    def term_annotate(self, args):
        return args[0]

    # Command Handlers
    def cmd_assert(self, args):
        self._assertions.append(args[0])

    def cmd_declare_const(self, args):
        name, sort = args
        sym = Symbol(name)
        sort = self._resolve_sort(sort)

        if sort == 'Real':
            self._assertions.append(Q.real(sym))
        elif sort == 'Int':
            self._assertions.append(Q.integer(sym))

        self._symbols[name] = sym

    def cmd_declare_fun(self, args):
        name = args[0]
        domain_sorts = args[1:-1]
        range_sort = args[-1]

        if not domain_sorts:
            # It's actually a constant
            return self.cmd_declare_const([name, range_sort])

        self._functions[name] = Function(name)

    def cmd_declare_sort(self, args):
        sort_name, arity = args
        self._sorts[sort_name] = int(arity)

    def cmd_declare_datatype(self, args):
        raise NotImplementedError(
            "Algebraic datatypes (declare-datatype) are not supported")

    def cmd_declare_datatypes(self, args):
        raise NotImplementedError(
            "Algebraic datatypes (declare-datatypes) are not supported")

    def cmd_define_fun(self, args):
        name = args[0]
        sorted_vars = args[1:-2]
        # range_sort = args[-2]
        body = args[-1]

        if not sorted_vars:
            self._symbols[name] = body
        else:
            var_names = [v[0] for v in sorted_vars]
            var_syms = [Symbol(v) for v in var_names]
            self._macros[name] = (var_syms, body)

    def cmd_define_fun_rec(self, args):
        raise NotImplementedError("Recursive macros are not supported in SymPy because AST generation evaluates bottom-up.")

    def cmd_define_const(self, args):
        name, _sort, body = args
        self._symbols[name] = body

    def cmd_define_sort(self, args):
        name = args[0]
        params = args[1:-1]
        sort = args[-1]
        if params:
            raise NotImplementedError(
                "Parametric sort definitions (define-sort with parameters) "
                "are not supported")
        self._sort_aliases[name] = self._resolve_sort(sort)

    def cmd_generic(self, args):
        raise UnknownSMTLibCommandError(f"Unknown SMT-LIB command: {args[0]}")

    def cmd_pop(self, args):
        levels = int(args[0]) if args else 1
        for _ in range(levels):
            if not self._assertion_stack:
                break
            self._assertions = self._assertion_stack.pop()
            self._symbols = self._symbol_stack.pop()
            self._macros = self._macro_stack.pop()
            self._functions = self._function_stack.pop()
            self._sorts = self._sort_stack.pop()
            self._sort_aliases = self._sort_alias_stack.pop()

    def cmd_push(self, args):
        levels = int(args[0]) if args else 1
        for _ in range(levels):
            self._assertion_stack.append(list(self._assertions))
            self._symbol_stack.append(dict(self._symbols))
            self._macro_stack.append(dict(self._macros))
            self._function_stack.append(dict(self._functions))
            self._sort_stack.append(dict(self._sorts))
            self._sort_alias_stack.append(dict(self._sort_aliases))

    def cmd_reset(self, args):
        self.__init__()

    def cmd_reset_assertions(self, args):
        self._assertions = []

    def cmd_set_logic(self, args):
        if args:
            self._logic = str(args[0])

    def cmd_set_info(self, args):
        if len(args) >= 2:
            self._info[str(args[0])] = args[1]

    def cmd_set_option(self, args):
        if len(args) >= 2:
            self._options[str(args[0])] = args[1]
