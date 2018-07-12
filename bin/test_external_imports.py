#!/usr/bin/env python
"""
Test that

from sympy import *

Doesn't import anything other than SymPy, it's hard dependencies (mpmath), and
hard optional dependencies (gmpy2, fastcache). Importing unnecessary libraries
can accidentally add hard dependencies to SymPy in the worst case, or at best
slow down the SymPy import time when they are installed.

Note, for this test to be effective, every external library that could
potentially be imported by SymPy must be installed.

TODO: Monkeypatch the importer to detect non-standard library imports even
when they aren't installed.

Based on code from
https://stackoverflow.com/questions/22195382/how-to-check-if-a-module-library-package-is-part-of-the-python-standard-library.
"""

# These libraries will always be imported with SymPy
hard_dependencies = ['mpmath']

# These libraries are optional, but are always imported at SymPy import time
# when they are installed. External libraries should only be added to this
# list if they are required for core SymPy functionality.
hard_optional_dependencies = ['gmpy', 'gmpy2', 'fastcache']

import sys
import os

stdlib = {p for p in sys.path if p.startswith(sys.prefix) and
    'site-packages' not in p}

existing_modules = list(sys.modules.keys())


# hook in-tree SymPy into Python path, if possible

this_path = os.path.abspath(__file__)
this_dir = os.path.dirname(this_path)
sympy_top = os.path.split(this_dir)[0]
sympy_dir = os.path.join(sympy_top, 'sympy')

if os.path.isdir(sympy_dir):
    sys.path.insert(0, sympy_top)

def test_external_imports():
    exec("from sympy import *", {})

    bad = []
    for mod in sys.modules:
        if '.' in mod and mod.split('.')[0] in sys.modules:
            # Only worry about the top-level modules
            continue

        if mod in existing_modules:
            continue

        if any(mod == i or mod.startswith(i + '.') for i in ['sympy'] +
            hard_dependencies + hard_optional_dependencies):
            continue

        if mod in sys.builtin_module_names:
            continue

        # if not hasattr(mod, '__file__'):
        #     # bad.append(mod)
        #     continue

        try:
            fname = sys.modules[mod].__file__
        except AttributeError:
            bad.append(mod)
            continue

        if fname.endswith(('__init__.py', '__init__.pyc', '__init__.pyo')):
            fname = os.path.dirname(fname)

        if os.path.dirname(fname) in stdlib:
            continue

        bad.append(mod)

    if bad:
        raise RuntimeError("""Unexpected external modules found when running 'from sympy import *':
    """ + '\n    '.join(bad))

    print("No unexpected external modules were imported with 'from sympy import *'!")

if __name__ == '__main__':
    test_external_imports()
