# Tests that require installed backends go into
# sympy/test_external/test_autowrap

import os
import tempfile
import shutil
from StringIO import StringIO

from sympy.utilities.autowrap import autowrap, binary_function, CythonCodeWrapper
from sympy.utilities.codegen import Routine, CCodeGen
from sympy.core import symbols

def get_string(dump_fn, routines, prefix="file", header=False, empty=False):
    """Wrapper for dump_fn. dump_fn writes its results to a stream object and
       this wrapper returns the contents of that stream as a string. This
       auxiliary function is used by many tests below.

       The header and the empty lines are not generator to facilitate the
       testing of the output.
    """
    output = StringIO()
    dump_fn(routines, output, prefix, header, empty)
    source = output.getvalue()
    output.close()
    return source

def test_cython_wrapper_scalar_function():
    x,y,z = symbols('xyz')
    expr = (x+y)*z
    routine = Routine("test", expr)
    code_gen = CythonCodeWrapper(CCodeGen())
    source = get_string(code_gen.dump_pyx, [routine])
    expected = (
            'cdef extern from "file.h":\n'
            '   double test(double x, double y, double z)\n'
            'def test_c(double x, double y, double z):\n'
            '   return test(x, y, z)\n'
    )
    assert source == expected

def test_cython_wrapper_outarg():
    from sympy import Equality
    x,y,z = symbols('xyz')
    code_gen = CythonCodeWrapper(CCodeGen())

    routine = Routine("test", Equality(z, x + y))
    source = get_string(code_gen.dump_pyx, [routine])
    expected = (
            'cdef extern from "file.h":\n'
            '   void test(double x, double y, double &z)\n'
            'def test_c(double x, double y):\n'
            '   cdef double z\n'
            '   test(x, y, z)\n'
            '   return z\n'
    )
    assert source == expected

def test_cython_wrapper_inoutarg():
    from sympy import Equality
    x,y,z = symbols('xyz')
    code_gen = CythonCodeWrapper(CCodeGen())
    routine = Routine("test", Equality(z, x + y + z))
    source = get_string(code_gen.dump_pyx, [routine])
    expected = (
            'cdef extern from "file.h":\n'
            '   void test(double x, double y, double &z)\n'
            'def test_c(double x, double y, double z):\n'
            '   test(x, y, z)\n'
            '   return z\n'
    )
    assert source == expected

def test_autowrap():
    x, y = symbols('x y')
    f = autowrap(x + y, backend='dummy')
    assert f() == str(x + y)

def test_autowrap_store_files():
    x, y = symbols('x y')
    tmp = tempfile.mkdtemp()
    try:
        f = autowrap(x + y, backend='dummy', tempdir=tmp)
        assert f() == str(x + y)
        assert os.access(tmp, os.F_OK)
    finally:
        shutil.rmtree(tmp)

def test_binary_function():
    x, y = symbols('x y')
    f = binary_function('f', x + y, backend='dummy')
    assert f._imp_() == str(x + y)
