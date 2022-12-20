import typing

import sympy
from sympy.core import Add, Mul
from sympy.core import Symbol, Expr, Float, Rational, Integer, Basic
from sympy.core.assumptions import assumptions
from sympy.core.function import UndefinedFunction, Function
from sympy.core.relational import Relational, Unequality, Equality, LessThan, GreaterThan, StrictLessThan, StrictGreaterThan
from sympy.functions.elementary.complexes import Abs
from sympy.functions.elementary.exponential import exp, log, Pow
from sympy.functions.elementary.hyperbolic import sinh, cosh, tanh
from sympy.functions.elementary.miscellaneous import Min, Max
from sympy.functions.elementary.piecewise import Piecewise
from sympy.functions.elementary.trigonometric import sin, cos, tan, asin, acos, atan, atan2
from sympy.logic.boolalg import And, Or, Xor, Implies, Boolean
from sympy.logic.boolalg import BooleanTrue, BooleanFalse, BooleanFunction, Not, ITE
from sympy.printing.printer import Printer
from sympy.sets import Interval


class SMTLibPrinter(Printer):
    printmethod = "_smtlib"

    # based on dReal, an automated reasoning tool for solving problems that can be encoded as first-order logic formulas over the real numbers.
    # dReal's special strength is in handling problems that involve a wide range of nonlinear real functions.
    _default_settings: dict = {
        'precision': None,
        'known_types': {
            bool: 'Bool',
            int: 'Int',
            float: 'Real'
        },
        'known_constants': {
            # pi: 'MY_VARIABLE_PI_DECLARED_ELSEWHERE',
        },
        'known_functions': {
            Add: '+',
            Mul: '*',

            Equality: '=',
            LessThan: '<=',
            GreaterThan: '>=',
            StrictLessThan: '<',
            StrictGreaterThan: '>',

            exp: 'exp',
            log: 'log',
            Abs: 'abs',
            sin: 'sin',
            cos: 'cos',
            tan: 'tan',
            asin: 'arcsin',
            acos: 'arccos',
            atan: 'arctan',
            atan2: 'arctan2',
            sinh: 'sinh',
            cosh: 'cosh',
            tanh: 'tanh',
            Min: 'min',
            Max: 'max',
            Pow: 'pow',

            And: 'and',
            Or: 'or',
            Xor: 'xor',
            Not: 'not',
            ITE: 'ite',
            Implies: '=>',
        }
    }

    symbol_table: dict

    def __init__(self, settings: dict = None, symbol_table=None):
        settings = settings or {}
        self.symbol_table = symbol_table or {}
        Printer.__init__(self, settings)
        self._precision = self._settings['precision']
        self._known_types = dict(self._settings['known_types'])
        self._known_constants = dict(self._settings['known_constants'])
        self._known_functions = dict(self._settings['known_functions'])

        for _ in self._known_types.values(): self._check_is_legal_name(_)
        for _ in self._known_constants.values(): self._check_is_legal_name(_)
        # for _ in self._known_functions.values(): self._check_is_legal_name(_)  # +, *, <, >, etc.

    def _is_legal_name(self, s: str):
        if not s: return False
        if s[0].isnumeric(): return False
        return all(_.isalnum() or _ == '_' for _ in s)

    def _check_is_legal_name(self, s: str):
        assert self._is_legal_name(s), f"Name `{s}` may not be legal in SMT-Lib."

    def _s_expr(self, op: str, args: typing.Union[list, tuple]) -> str:
        args_str = ' '.join(
            a if isinstance(a, str)
            else self._print(a)
            for a in args
        )
        return f'({op} {args_str})'

    def _print_Function(self, e):
        if e in self._known_functions:
            op = self._known_functions[e]
        elif type(e) in self._known_functions:
            op = self._known_functions[type(e)]
        elif type(type(e)) == UndefinedFunction:
            op = e.name
        else:
            op = self._known_functions[e]  # throw KeyError

        return self._s_expr(op, e.args)

    def _print_Relational(self, e: Relational):
        return self._print_Function(e)

    def _print_BooleanFunction(self, e: BooleanFunction):
        return self._print_Function(e)

    def _print_Expr(self, e: Expr):
        return self._print_Function(e)

    def _print_Unequality(self, e: Unequality):
        if type(e) in self._known_functions:
            return self._print_Relational(e)  # default
        else:
            eq_op = self._known_functions[Equality]
            not_op = self._known_functions[Not]
            return self._s_expr(not_op, [self._s_expr(eq_op, e.args)])

    def _print_Piecewise(self, e: Piecewise):
        def _print_Piecewise_recursive(args: typing.Union[list, tuple]):
            e, c = args[0]
            if len(args) == 1:
                assert (c is True) or isinstance(c, BooleanTrue), "Piecewise expression must end in (expr, True) statement."
                return self._print(e)
            else:
                ite = self._known_functions[ITE]
                return self._s_expr(ite, [
                    c, e, _print_Piecewise_recursive(args[1:])
                ])

        return _print_Piecewise_recursive(e.args)

    def _print_Interval(self, e: Interval):
        if e.start.is_infinite and e.end.is_infinite:
            return ''
        elif e.start.is_infinite != e.end.is_infinite:
            raise ValueError(f'One-sided intervals (`{e}`) are not supported in SMT.')
        else:
            return f'[{e.start}, {e.end}]'

    # todo: Sympy does not support quantifiers yet as of 2022, but quantifiers can be handy in SMT.
    # For now, users can extend this class and build in their own quantifier support.
    # See `test_quantifier_extensions()` in test_smtlib.py for an example of how this might look.

    # def _print_ForAll(self, e: ForAll):
    #     return self._s('forall', [
    #         self._s('', [
    #             self._s(sym.name, [self._type_name(sym), Interval(start, end)])
    #             for sym, start, end in e.limits
    #         ]),
    #         e.function
    #     ])

    def _print_BooleanTrue(self, x: BooleanTrue):
        return 'true'

    def _print_BooleanFalse(self, x: BooleanFalse):
        return 'false'

    def _print_Float(self, x: Float):
        f = x.evalf(self._precision) if self._precision else x.evalf()
        return str(f).rstrip('0')

    def _print_float(self, x: float):
        return str(x)

    def _print_Rational(self, x: Rational):
        return self._s_expr('/', [x.p, x.q])

    def _print_Integer(self, x: Integer):
        assert x.q == 1
        return str(x.p)

    def _print_int(self, x: int):
        return str(x)

    def _print_Symbol(self, x: Symbol):
        self._check_is_legal_name(x.name)
        return x.name

    def _print_RandomSymbol(self, x):
        self._check_is_legal_name(x.name)
        return x.name

    def _print_NumberSymbol(self, x):
        name = self._known_constants.get(x)
        return name if name else self._print_Float(x)

    def _print_UndefinedFunction(self, x):
        self._check_is_legal_name(x.name)
        return x.name

    def _print_Exp1(self, x):
        return (
            self._print_Function(exp(1, evaluate=False))
            if exp in self._known_functions else
            self._print_NumberSymbol(x)
        )

    def emptyPrinter(self, expr):
        raise NotImplementedError(f'Cannot convert `{repr(expr)}` of type `{type(expr)}` to SMT.')


def smtlib_code(
    expr,
    auto_assert=True, auto_declare=True,
    precision=None,
    symbol_table=None,
    known_types=None, known_constants=None, known_functions=None,
    prefix_expressions=None, suffix_expressions=None,
    log_warn=None
):
    r"""Converts ``expr`` to a string of smtlib code.

    Parameters
    ==========

    expr : Expr | List[Expr]
        A SymPy expression or system to be converted.
    auto_assert : bool, optional
        If false, do not modify expr and produce only the S-Expression equivalent of expr.
        If true, assume expr is a system and assert each boolean element.
    auto_declare : bool, optional
        If false, do not produce declarations for the symbols used in expr.
        If true, prepend all necessary declarations for variables used in expr based on symbol_table.
    precision : integer, optional
        The ``evalf(..)`` precision for numbers such as pi.
    symbol_table : dict, optional
        A dictionary where keys are ``Symbol`` or ``Function`` instances and values are their Python type i.e. ``bool``, ``int``, ``float``, or ``Callable[...]``.
        If incomplete, an attempt will be made to infer types from ``expr``.
    known_types: dict, optional
        A dictionary where keys are ``bool``, ``int``, ``float`` etc. and values are their corresponding SMT type names.
        If not given, a partial listing compatible with several solvers will be used.
    known_functions : dict, optional
        A dictionary where keys are ``Function``, ``Relational``, ``BooleanFunction``, or ``Expr`` instances and values are their SMT string representations.
        If not given, a partial listing optimized for dReal solver (but compatible with others) will be used.
    known_constants: dict, optional
        A dictionary where keys are ``NumberSymbol`` instances and values are their SMT variable names.
        When using this feature, extra caution must be taken to avoid naming collisions between user symbols and listed constants.
        If not given, constants will be expanded inline i.e. ``3.14159`` instead of ``MY_SMT_VARIABLE_FOR_PI``.
    prefix_expressions: list, optional
        A list of lists of ``str`` and/or expressions to convert into SMTLib and prefix to the output.
    suffix_expressions: list, optional
        A list of lists of ``str`` and/or expressions to convert into SMTLib and postfix to the output.
    log_warn: lambda function, optional
        A function to record all warnings during potentially risky operations.
        Soundness is a core value in SMT solving, so it is good to log all assumptions made.

    Examples
    ========
    >>> from sympy import smtlib_code, symbols, sin, Eq
    >>> x = symbols('x')
    >>> smtlib_code(sin(x).series(x).removeO(), log_warn=print)
    Could not infer type of `x`. Defaulting to float.
    Non-Boolean expression `x**5/120 - x**3/6 + x` will not be asserted. Converting to SMTLib verbatim.
    '(declare-const x Real)\n(+ x (* (/ -1 6) (pow x 3)) (* (/ 1 120) (pow x 5)))'

    >>> from sympy import Rational
    >>> x, y, tau = symbols("x, y, tau")
    >>> smtlib_code((2*tau)**Rational(7, 2), log_warn=print)
    Could not infer type of `tau`. Defaulting to float.
    Non-Boolean expression `8*sqrt(2)*tau**(7/2)` will not be asserted. Converting to SMTLib verbatim.
    '(declare-const tau Real)\n(* 8 (pow 2 (/ 1 2)) (pow tau (/ 7 2)))'

    ``Piecewise`` expressions are implemented with ``ite`` expressions by default.
    Note that if the ``Piecewise`` lacks a default term, represented by
    ``(expr, True)`` then an error will be thrown.  This is to prevent
    generating an expression that may not evaluate to anything.

    >>> from sympy import Piecewise
    >>> pw = Piecewise((x + 1, x > 0), (x, True))
    >>> smtlib_code(Eq(pw, 3), symbol_table={x: float}, log_warn=print)
    '(declare-const x Real)\n(assert (= (ite (> x 0) (+ 1 x) x) 3))'

    Custom printing can be defined for certain types by passing a dictionary of
    PythonType : "SMT Name" to the ``known_types``, ``known_constants``, and ``known_functions`` kwargs.

    >>> from typing import Callable
    >>> from sympy import Function, Add
    >>> f = Function('f')
    >>> g = Function('g')
    >>> smt_builtin_funcs = {  # functions our SMT solver will understand
    ...   f: "existing_smtlib_fcn",
    ...   Add: "sum",
    ... }
    >>> user_def_funcs = {  # functions defined by the user must have their types specified explicitly
    ...   g: Callable[[int], float],
    ... }
    >>> smtlib_code(f(x) + g(x), symbol_table=user_def_funcs, known_functions=smt_builtin_funcs, log_warn=print)
    Non-Boolean expression `f(x) + g(x)` will not be asserted. Converting to SMTLib verbatim.
    '(declare-const x Int)\n(declare-fun g (Int) Real)\n(sum (existing_smtlib_fcn x) (g x))'
    """
    log_warn = log_warn or (lambda _: None)

    if not isinstance(expr, list): expr = [expr]
    expr = [
        sympy.sympify(_, strict=True, evaluate=False, convert_xor=False)
        for _ in expr
    ]

    if not symbol_table: symbol_table = {}
    symbol_table = _auto_infer_smtlib_types(
        expr, symbol_table, auto_assert
    )
    # See [FALLBACK RULES]
    # Need SMTLibPrinter to populate known_functions and known_constants first.

    settings = {}
    if precision: settings['precision'] = precision
    del precision

    if known_types: settings['known_types'] = known_types
    del known_types

    if known_functions: settings['known_functions'] = known_functions
    del known_functions

    if known_constants: settings['known_constants'] = known_constants
    del known_constants

    if not prefix_expressions: prefix_expressions = []
    if not suffix_expressions: suffix_expressions = []

    p = SMTLibPrinter(settings, symbol_table)
    del symbol_table

    # [FALLBACK RULES]
    for e in expr:
        for sym in _atoms_symbols_preserve_rv(e):
            if (
                sym.is_symbol and
                sym not in p._known_constants and
                sym not in p.symbol_table
            ):
                log_warn(f"Could not infer type of `{sym}`. Defaulting to float.")
                p.symbol_table[sym] = float
        for fun in e.atoms(Function):
            if (
                fun.is_Function and
                type(fun) not in p._known_functions and
                type(fun) not in p.symbol_table and
                not fun.is_Piecewise
            ): raise TypeError(
                f"Unknown type of undefined function `{fun}`. "
                f"Must be mapped to ``str`` in known_functions or mapped to ``Callable[..]`` in symbol_table."
            )

    declarations = []
    if auto_declare:
        declarations = {
                           sym.name: sym for e in expr for sym in e.free_symbols  # .free_symbols preserves random variables
                           if sym not in p._known_constants
                       } | {
                           fnc.name: fnc for e in expr for fnc in e.atoms(Function)
                           if type(fnc) not in p._known_functions and not fnc.is_Piecewise
                       }
        declarations = [
            _auto_declare_smtlib(sym_or_func, p, log_warn) for sym_or_func in declarations.values()
        ]

    if auto_assert:
        expr = [_auto_assert_smtlib(e, p, log_warn) for e in expr]

    # return SMTLibPrinter().doprint(expr)
    return '\n'.join([
        # ';; PREFIX EXPRESSIONS',
        *[
            e if isinstance(e, str) else p.doprint(e)
            for e in prefix_expressions
        ],

        # ';; DECLARATIONS',
        *[line for block in sorted([_ for _ in declarations if _], key=lambda b: b[0]) for line in block],

        # ';; EXPRESSIONS',
        *[
            e if isinstance(e, str) else p.doprint(e)
            for e in expr
        ],

        # ';; SUFFIX EXPRESSIONS',
        *[
            e if isinstance(e, str) else p.doprint(e)
            for e in suffix_expressions
        ],
    ])


def _atoms_symbols_preserve_rv(expr):
    from sympy.stats.rv import RandomSymbol
    random_symbols = expr.atoms(RandomSymbol)
    simple_symbols = set(expr.atoms(Symbol)) - {_.symbol for _ in random_symbols}
    return list(simple_symbols) + list(random_symbols)


def _auto_declare_smtlib(sym: typing.Union[Symbol, Function], p: SMTLibPrinter, log_warn: typing.Callable[[str], None]):
    if sym.is_symbol:
        type_signature = p.symbol_table[sym]
        assert isinstance(type_signature, type)
        type_signature = p._known_types[type_signature]

        from sympy.stats.rv import RandomSymbol
        current_assumptions = assumptions(sym) | {'__is_random_symbol': isinstance(sym, RandomSymbol)}
        unsupported_assumptions = {
            'infinite', 'antihermitian', 'transcendental', 'imaginary', 'irrational', 'noninteger', 'even', 'odd', 'prime', 'composite'
        }
        supported_assumptions = {
            'zero': lambda: Equality(sym, 0, evaluate=False),
            'nonzero': lambda: Unequality(sym, 0, evaluate=False),
            'positive': lambda: StrictGreaterThan(sym, 0, evaluate=False),
            'nonnegative': lambda: GreaterThan(sym, 0, evaluate=False),
            'negative': lambda: StrictLessThan(sym, 0, evaluate=False),
            'nonpositive': lambda: LessThan(sym, 0, evaluate=False),
            '__is_random_symbol': lambda: sym.pspace.domain.as_boolean().simplify()
        }
        for a in unsupported_assumptions:
            if current_assumptions.get(a) and not current_assumptions.get('zero'):  # zero checks pretty much everything
                raise ValueError(f"Cannot automatically assert '{a}'-ness when declaring `{sym}`. Please assert explicitly.")
        return [
                   p._s_expr('declare-const', [sym, type_signature])
               ] + [
                   p._s_expr('assert', [predicate()])
                   for a, predicate in supported_assumptions.items()
                   if current_assumptions.get(a)
               ]

    elif sym.is_Function:
        type_signature = p.symbol_table[type(sym)]
        assert callable(type_signature)
        type_signature = [p._known_types[_] for _ in type_signature.__args__]
        assert len(type_signature) > 0
        params_signature = f"({' '.join(type_signature[:-1])})"
        return_signature = type_signature[-1]
        return [p._s_expr('declare-fun', [type(sym), params_signature, return_signature])]

    else:
        raise ValueError(f"Non-Symbol/Function `{sym}` will not be declared.")


def _auto_assert_smtlib(e: Expr, p: SMTLibPrinter, log_warn: typing.Callable[[str], None]):
    if isinstance(e, Boolean) or (
        e in p.symbol_table and p.symbol_table[e] == bool
    ) or (
        e.is_Function and
        type(e) in p.symbol_table and
        p.symbol_table[type(e)].__args__[-1] == bool
    ):
        return p._s_expr('assert', [e])
    else:
        log_warn(f"Non-Boolean expression `{e}` will not be asserted. Converting to SMTLib verbatim.")
        return e


def _auto_infer_smtlib_types(
    exprs: typing.List[Basic],
    symbol_table: dict = None,
    auto_assert: bool = True
) -> dict:
    # [TYPE INFERENCE RULES]
    # X is alone in an expr and auto-assert => X is bool
    # X in BooleanFunction.args => X is bool
    # X matches to a bool param of a symbol_table function => X is bool
    # X matches to an int param of a symbol_table function => X is int
    # X.is_integer => X is int
    # X is random and X.pspace is continuous => X is float
    # X == Y, where X is T => Y is T

    # [FALLBACK RULES]
    # see _auto_declare_smtlib(..)
    # X is not bool and X is not int and X is Symbol => X is float
    # else (e.g. X is Function) => error. must be specified explicitly.

    _symbols = dict(symbol_table) if symbol_table else {}

    def safe_update(syms: set, inf):
        for s in syms:
            assert s.is_symbol
            if (old_type := _symbols.setdefault(s, inf)) != inf:
                raise TypeError(f"Could not infer type of `{s}`. Apparently both `{old_type}` and `{inf}`?")

    # EXPLICIT TYPES
    safe_update({
        e
        for e in exprs
        if e.is_symbol and auto_assert
    }, bool)

    safe_update({
        symbol
        for e in exprs
        for boolfunc in e.atoms(BooleanFunction)
        for symbol in boolfunc.args
        if symbol.is_symbol
    }, bool)

    safe_update({
        symbol
        for e in exprs
        for boolfunc in e.atoms(Function)
        if type(boolfunc) in _symbols
        for symbol, param in zip(boolfunc.args, _symbols[type(boolfunc)].__args__)
        if symbol.is_symbol and param == bool
    }, bool)

    safe_update({
        symbol
        for e in exprs
        for intfunc in e.atoms(Function)
        if type(intfunc) in _symbols
        for symbol, param in zip(intfunc.args, _symbols[type(intfunc)].__args__)
        if symbol.is_symbol and param == int
    }, int)

    safe_update({
        symbol
        for e in exprs
        for symbol in _atoms_symbols_preserve_rv(e)
        if symbol.is_integer
    }, int)

    safe_update({
        symbol
        for e in exprs
        for symbol in _atoms_symbols_preserve_rv(e)
        if symbol.is_real and not symbol.is_integer
    }, float)

    # CONTINUOUS RANDOM VARIABLE RULE
    from sympy.stats.rv import RandomSymbol
    safe_update({
        rv
        for e in exprs
        for rv in e.atoms(RandomSymbol)
        if rv.pspace.is_Continuous
    }, float)

    # EQUALITY RELATION RULE
    rels = [rel for expr in exprs for rel in expr.atoms(Equality)]
    rels = [
               (rel.lhs, rel.rhs) for rel in rels if rel.lhs.is_symbol
           ] + [
               (rel.rhs, rel.lhs) for rel in rels if rel.rhs.is_symbol
           ]
    for infer, reltd in rels:
        inference = (
            _symbols[infer] if infer in _symbols else
            _symbols[reltd] if reltd in _symbols else

            _symbols[type(reltd)].__args__[-1]
            if reltd.is_Function and type(reltd) in _symbols else

            bool if reltd.is_Boolean else
            int if reltd.is_integer or reltd.is_Integer else
            float if reltd.is_real else
            None
        )
        if inference: safe_update({infer}, inference)

    return _symbols
