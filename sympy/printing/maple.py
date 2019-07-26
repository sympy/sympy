"""
Maple code printer

The MapleCodePrinter converts single sympy expressions into single
Maple expressions, using the functions defined in the Maple objects where possible.


FIXME: This module is still under actively developed. Some functions may be not completed.
"""

from __future__ import print_function, division

from sympy.codegen.ast import Assignment
from sympy.core import S
from sympy.core.compatibility import string_types, range
from sympy.printing.codeprinter import CodePrinter
from sympy.printing.precedence import precedence, PRECEDENCE

_known_func_same_name = [
    'sin', 'cos', 'tan', 'sec', 'csc', 'cot', 'sinh', 'cosh', 'tanh', 'sech',
    'csch', 'coth', 'exp', 'floor'
]

known_functions = {
    # Sympy -> Maple
    # FIXME: maybe need added
    'Abs': 'abs',
    'log': 'ln',
    'asin': 'arcsin',
    'acos': 'arccos',
    'atan': 'arctan',
    'asec': 'arcsec',
    'acsc': 'arccsc',
    'acot': 'arccot',
    'asinh': 'arcsinh',
    'acosh': 'arccosh',
    'atanh': 'arctanh',
    'asech': 'arcsech',
    'acsch': 'arccsch',
    'acoth': 'arccoth',
    'ceiling': 'ceil',
}

atomic_expr = {
    # Sympy -> Maple
    S('pi'): 'Pi',
    S('I'): 'I',
    S('oo'): 'infinity'
}

for _func in _known_func_same_name:
    known_functions[_func] = _func


class MapleCodePrinter(CodePrinter):
    """
    Printer which converts single sympy expressions into single.
    """
    def __init__(self, settings=None):
        if settings is None:
            settings = dict()
        CodePrinter.__init__(settings)
        self.known_functions = dict(known_functions)
        userfuncs = settings.get('user_functions', {})
        self.known_functions.update(userfuncs)

    def _get_statement(self, codestring):
        return "%s;" % codestring

    def _get_comment(self, text):
        return "# {0}".format(text)

    def _declare_number_const(self, name, value):
        return "{0} := {1};".format(name,
                                    value.evalf(self._settings['precision']))

    def _print_Pow(self, expr):
        PREC = precedence(expr)
        if expr.exp == -1:
            return '1/%s' % (self.parenthesize(expr.base, PREC))
        elif expr.exp == 0.5 or expr.exp == S(1) / 2:
            return 'sqrt(%s)' % self._print(expr.base)
        elif expr.exp == -0.5 or expr.exp == -S(1) / 2:
            return '1/sqrt(%s)' % self._print(expr.base)
        else:
            return '%s^%s' % (self._print(expr.base), self._print(expr.exp))

    def _print_Rational(self, expr):
        p, q = int(expr.p), int(expr.q)
        return '%d/%d' % (p, q)

    def _print_Relational(self, expr):
        lhs_code = self._print(expr.lhs)
        rhs_code = self._print(expr.rhs)
        op = expr.rel_op
        # FIXME: rel_op might be different from maple.
        return "{0} {1} {2}".format(lhs_code, op, rhs_code)

    def _print_AtomicExpr(self, expr):
        """
        Print constant like pi, e ...
        """
        return atomic_expr[expr]

    def _print_Idx(self, expr):
        return self._print(expr.label)

    def _print_MatrixBase(self, expr):
        if expr.cols == 0 or expr.rows == 0:
            return '<>'
        elif expr.cols == 1:
            return '<%s>' % ','.join(list(expr.col(0)))
        elif expr.rows == 1:
            return '<%s>' % '|'.join(list(expr.col(0)))
        else:
            _row_content_list = [
                '<%s>' % '|'.join(list(expr.row(i))) for i in range(expr.rows)
            ]
            return '<%s>' % ','.join(_row_content_list)


def maple_code(expr, assign_to=None, **settings):
    r"""Converts `expr` to a string of Julia code.

        Parameters
        ==========

        expr : Expr
            A sympy expression to be converted.
        assign_to : optional
            When given, the argument is used as the name of the variable to which
            the expression is assigned.  Can be a string, ``Symbol``,
            ``MatrixSymbol``, or ``Indexed`` type.  This can be helpful for
            expressions that generate multi-line statements.
        precision : integer, optional
            The precision for numbers such as pi  [default=16].
        user_functions : dict, optional
            A dictionary where keys are ``FunctionClass`` instances and values are
            their string representations.  Alternatively, the dictionary value can
            be a list of tuples i.e. [(argument_test, cfunction_string)].  See
            below for examples.
        human : bool, optional
            If True, the result is a single string that may contain some constant
            declarations for the number symbols.  If False, the same information is
            returned in a tuple of (symbols_to_declare, not_supported_functions,
            code_text).  [default=True].
        contract: bool, optional
            If True, ``Indexed`` instances are assumed to obey tensor contraction
            rules and the corresponding nested loops over indices are generated.
            Setting contract=False will not generate loops, instead the user is
            responsible to provide values for the indices in the code.
            [default=True].
        inline: bool, optional
            If True, we try to create single-statement code instead of multiple
            statements.  [default=True].

        FIXME: examples should be added.
    """
    return MapleCodePrinter(settings).doprint(expr, assign_to)


def print_maple_code(expr, **settings):
    """Prints the Maple representation of the given expression.

        See `maple_code` for the meaning of the optional arguments.
    """
    print(maple_code(expr, settings))
