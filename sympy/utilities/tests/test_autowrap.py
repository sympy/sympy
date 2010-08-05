# Tests that require installed backends go into
# sympy/test_external/test_autowrap

import os
import tempfile
import shutil

from sympy.utilities.autowrap import autowrap, binary_function
from sympy.core import symbols

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
