from __future__ import annotations

from sympy import (
    Symbol, Function, Integer, Float, Eq, Ne, Not, And, Or,
    Implies, Xor, Piecewise, Add, Mul, StrictLessThan, StrictGreaterThan,
    LessThan, GreaterThan, Expr
)
from sympy.assumptions import Q
from sympy.functions.elementary.complexes import Abs
from sympy.external import import_module
from sympy.parsing.smtlib.lark.smtlib_parser import UnknownSMTLibCommandError

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
        self._options = {}
        self._info = {}

        self._extended_ops = {
            'to_real', 'to_int', 'is_int',
            'concat', 'bvnot', 'bvand', 'bvor', 'bvneg', 'bvadd', 'bvmul',
            'bvudiv', 'bvurem', 'bvshl', 'bvlshr', 'bvsub', 'bvult', 'bvxor',
            'bvnand', 'bvnor', 'bvxnor', 'bvcomp', 'bvsdiv', 'bvsrem',
            'bvsmod', 'bvashr', 'bvule', 'bvugt', 'bvuge', 'bvslt', 'bvsle',
            'bvsgt', 'bvsge', 'bv2nat', 'bv2int', 'ubv_to_int',
            'ext_rotate_left', 'ext_rotate_right', 'str.len', 'str.++',
            'str.at', 'str.contains', 'str.indexof', 'str.replace',
            'str.substr', 'str.prefixof', 'str.suffixof', 'str.to.int',
            'int.to.str', 'select', 'store', 'div', 'mod'
        }

        self._assertion_stack = []
        self._symbol_stack = []
        self._macro_stack = []
        self._function_stack = []
        self._sort_stack = []

    def get_result(self):
        return self._symbols, self._assertions

    # Tokens
    def SYMBOL(self, token):
        return str(token)

    def QUOTED_SYMBOL(self, token):
        return str(token)

    def KEYWORD(self, token):
        return str(token)

    def NUMERAL(self, token):
        return Integer(token)

    def DECIMAL(self, token):
        return Float(token)

    def HEXADECIMAL(self, token):
        return Integer(int(str(token).replace('#x', '0x'), 0))

    def BINARY(self, token):
        return Integer(int(str(token).replace('#b', '0b'), 0))

    def STRING(self, token):
        # SMT-LIB strings are quoted. We return them as a Symbol because
        # SymPy doesn't have a native String type.
        return Symbol(str(token))

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

    def qual_identifier(self, args):
        # Can be just symbol
        if len(args) == 1:
            return args[0]
        # Or indexed: ("_", symbol, numerals...)
        if args[0] == '_':
            return Function("_")(args[1], *args[2:])
        # Or as: ("as", symbol, sort) or ("as", ("_", ...), sort)
        if args[0] == 'as':
            # SymPy doesn't explicitly track sort for expressions in AST,
            # so we just return the underlying identifier.
            if args[1] == '_':
                return Function("_")(args[2], *args[3:-1])
            return args[1]

    def term_qual_identifier(self, args):
        atom = args[0]
        # If atom is an indexed identifier (Function), return it directly
        if isinstance(atom, (Function, Expr)):
            return atom

        if atom in self._symbols:
            return self._symbols[atom]

        atom_str = str(atom).lower()
        if atom_str == 'true':
            return True
        if atom_str == 'false':
            return False

        return Symbol(str(atom))

    def term_apply(self, args):
        # args[0] is qual_identifier
        op_node = args[0]
        params = args[1:]

        # If the operator is an indexed identifier, pass it through directly
        # as a Function call
        if isinstance(op_node, (Function, Expr)):
            return op_node(*params)

        op = str(op_node)

        if op in self._macros:
            macro_args, macro_body = self._macros[op]
            return macro_body.subs(dict(zip(macro_args, params)))

        if op in self._functions:
            return self._functions[op](*params)

        op_lower = op.lower()

        if op_lower == 'not':
            return Not(params[0])
        elif op_lower == 'and':
            return And(*params)
        elif op_lower == 'or':
            return Or(*params)
        elif op_lower == '=>':
            return Implies(params[0], params[1])
        elif op_lower == 'xor':
            return Xor(*params)
        elif op_lower == '=':
            if len(params) == 2:
                return Eq(params[0], params[1])
            return And(*[Eq(params[i], params[i+1]) for i in range(len(params)-1)])
        elif op_lower == 'distinct':
            if len(params) == 2:
                return Ne(params[0], params[1])
            exprs = []
            for i in range(len(params)):
                for j in range(i + 1, len(params)):
                    exprs.append(Ne(params[i], params[j]))
            return And(*exprs)
        elif op_lower == '<':
            return StrictLessThan(params[0], params[1])
        elif op_lower == '<=':
            return LessThan(params[0], params[1])
        elif op_lower == '>':
            return StrictGreaterThan(params[0], params[1])
        elif op_lower == '>=':
            return GreaterThan(params[0], params[1])
        elif op_lower == '+':
            return Add(*params)
        elif op_lower == '-':
            if len(params) == 1:
                return -params[0]
            return params[0] - Add(*params[1:])
        elif op_lower == '*':
            return Mul(*params)
        elif op_lower == '/':
            if len(params) == 2:
                return params[0] / params[1]
            return params[0] / Mul(*params[1:])
        elif op_lower == 'ite':
            return Piecewise((params[1], params[0]), (params[2], True))
        elif op_lower == 'abs':
            return Abs(params[0])

        if op_lower in self._extended_ops:
            return Function(op)(*params)

        raise UnknownSMTLibCommandError(f"Unknown operator or function: {op}")


    def term_let(self, args):
        bindings = args[:-1]
        term = args[-1]
        for var, val in bindings:
            term = term.subs(Symbol(var), val)
        return term

    def term_forall(self, args):
        vars_list = args[:-1]
        body = args[-1]
        # vars_list comes from sorted_var nodes, which are (symbol, sort)
        var_syms = tuple(v[0] for v in vars_list)
        return Function('forall')(var_syms, body)

    def term_exists(self, args):
        vars_list = args[:-1]
        body = args[-1]
        var_syms = tuple(v[0] for v in vars_list)
        return Function('exists')(var_syms, body)

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
        # We handle this same as define-fun for now, SymPy will leave
        # unresolvable recursive invocations as uninterpreted.
        self.cmd_define_fun(args)

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

    def cmd_push(self, args):
        levels = int(args[0]) if args else 1
        for _ in range(levels):
            self._assertion_stack.append(list(self._assertions))
            self._symbol_stack.append(dict(self._symbols))
            self._macro_stack.append(dict(self._macros))
            self._function_stack.append(dict(self._functions))
            self._sort_stack.append(dict(self._sorts))

    def cmd_reset(self, args):
        self.__init__()

    def cmd_reset_assertions(self, args):
        self._assertions = []
