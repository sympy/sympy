# This file contains tests that exercise multiple AST nodes

from sympy.external import import_module
from sympy.printing.ccode import ccode
from sympy.utilities._compilation import compile_link_import_strings
from sympy.utilities._compilation.util import TemporaryDirectory
from sympy.utilities.pytest import skip
from sympy.sets import Range
from sympy.codegen.ast import (
    FunctionDefinition, FunctionPrototype, Variable, Pointer, real, Assignment,
    integer, Variable, CodeBlock, While, Statement
)
from sympy.codegen.cnodes import void, PreIncrement
from sympy.codegen.cutils import render_as_source_file

np = import_module('numpy')

def _mk_func1():
    declars = n, inp, out = Variable('n', integer), Pointer('inp', real), Pointer('out', real)
    i = Variable('i', integer)
    whl = While(i<n, [Assignment(out[i], inp[i]), Statement(PreIncrement(i))])
    body = CodeBlock(Statement(i.as_Declaration(value=0)), whl)
    return FunctionDefinition(void, 'our_test_function', declars, body)


def _render_compile_import(funcdef, build_dir):
    code_str = render_as_source_file(funcdef, settings=dict(contract=False))
    declar = ccode(FunctionPrototype.from_FunctionDefinition(funcdef))
    return compile_link_import_strings([
        ('our_test_func.c', code_str),
        ('_our_test_func.pyx', ("cdef extern {declar}\n"
                                "def _{fname}({typ}[:] inp, {typ}[:] out):\n"
                                "    {fname}(inp.size, &inp[0], &out[0])").format(
                                    declar=declar, fname=funcdef.name, typ='double'
                                ))
    ], build_dir=build_dir)


def test_copying_function():
    if not np:
        skip("numpy not installed.")

    info = None
    with TemporaryDirectory() as folder:
        mod, info = _render_compile_import(_mk_func1(), build_dir=folder)
        inp = np.arange(10.0)
        out = np.empty_like(inp)
        mod._our_test_function(inp, out)
        assert np.allclose(inp, out)
