from sympy.codegen.lfort import sympy_to_lfortran
from sympy.external import import_module
from sympy.printing.fcode import fcode

from sympy.utilities.pytest import XFAIL

from sympy import Integer, Float

lfortran = import_module("lfortran")


class ModifiedVisitor(lfortran.semantic.analyze.ExprVisitor):

    types = {
        'real': lfortran.semantic.analyze.Real(),
        'integer': lfortran.semantic.analyze.Integer()
    }

    def __init__(self, global_scope, expr_type):
        self.curr_type = self.types[expr_type]
        super().__init__(global_scope)

    # NOTE: We need to override this method because the analyzer currently
    # checks to see if a variable is bound within a scope. We'll be testing
    # against single expressions (such as "x + 1"), so that definitely won't be
    # the case
    def visit_Name(self, node):
        node._type = self.curr_type


# Adapted from LFortran's test_parser.py
def to_tuple(t):
    if t is None or isinstance(t, (str, int, complex, float)):
        return t
    elif isinstance(t, list):
        return [to_tuple(e) for e in t]
    result = [t.__class__.__name__]
    if t._fields is None:
        return tuple(result)
    elif isinstance(t, lfortran.ast.ast.Num):
        result.append(to_tuple(t.o))
    else:
        for f in t._fields:
            result.append(to_tuple(getattr(t, f)))
    return tuple(result)


def ast_to_tuple(ast, expr_type):
    """Parse a Fortran expression, producing a tuple from the annotated AST."""
    symbol_table = lfortran.semantic.analyze.create_symbol_table(ast)
    visitor = ModifiedVisitor(symbol_table, expr_type)
    visitor.visit(ast)
    return to_tuple(ast)


def lfort_tup(expr, expr_type):
    out_asr = sympy_to_lfortran(expr)
    return ast_to_tuple(lfortran.asr_to_ast(out_asr), expr_type)


def fcode_tup(expr, expr_type):
    return ast_to_tuple(lfortran.src_to_ast(fcode(expr),
                                            translation_unit=False), expr_type)


def test_integer():
    exprs = [
        Integer(1),
        Integer(0),
        Integer(42),
        Integer(1000000),
    ]
    expr_type = "integer"
    for expr in exprs:
        assert lfort_tup(expr, expr_type) == fcode_tup(expr, expr_type)


@XFAIL
def test_integer_failing():
    exprs = [
        Integer(-42),       # Parsed as unary subtraction, generated as BinOp
        Integer(-1)         # Node visitor not implemented
    ]
    expr_type = "integer"
    for expr in exprs:
        assert lfort_tup(expr, expr_type) == fcode_tup(expr, expr_type)


def test_float():
    exprs = [
        Float(0.0),
        Float(1e20),
    ]
    expr_type = "real"
    for expr in exprs:
        assert lfort_tup(expr, expr_type) == fcode_tup(expr, expr_type)


def test_symbol_addition():
    from sympy.abc import x, y, z
    exprs = [
        x + y,
        x + 1.0,
        1.0 * x + 3.0,
        x + y + z
    ]
    expr_type = "real"
    for expr in exprs:
        assert lfort_tup(expr, expr_type) == fcode_tup(expr, expr_type)
