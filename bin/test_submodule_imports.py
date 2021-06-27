#!/usr/bin/env python
"""
Test that

from sympy import *

only imports those sympy submodules that have names that are part of the
top-level namespace.
"""

import sys
import os

# hook in-tree SymPy into Python path, if possible

this_path = os.path.abspath(__file__)
this_dir = os.path.dirname(this_path)
sympy_top = os.path.split(this_dir)[0]
sympy_dir = os.path.join(sympy_top, 'sympy')

if os.path.isdir(sympy_dir):
    sys.path.insert(0, sympy_top)

submodule_whitelist = [
    'algebras',
    'assumptions',
    'calculus',
    'concrete',
    'core',
    'deprecated',
    'discrete',
    'external',
    'functions',
    'geometry',
    'integrals',
    'interactive',
    'logic',
    'matrices',
    'multipledispatch',
    'ntheory',
    'parsing',
    'plotting',
    'polys',
    'printing',
    'release',
    'series',
    'sets',
    'simplify',
    'solvers',
    'strategies',
    'tensor',
    'testing',
    'utilities',
]

def test_submodule_imports():
    if 'sympy' in sys.modules:
        raise RuntimeError("SymPy has already been imported, the test_submodule_imports test cannot run")

    exec("from sympy import *", {})

    for mod in sys.modules:
        if not mod.startswith('sympy'):
            continue

        if not mod.count('.') == 1:
            continue

        _, submodule = mod.split('.')

        if submodule not in submodule_whitelist:
            sys.exit(f"""\
Error: The submodule {mod} was imported with 'from sympy import *', but it was
not expected to be.

If {mod} is a new module that has functions that are imported at the
top-level, then the whitelist in bin/test_submodule_imports should be updated.
If it is not, the place that imports it should be modified so that it does not
get imported at the top-level, e.g., by moving the 'import {mod}' import
inside the function that uses it.

If you are unsure which code is importing {mod}, it may help to add 'raise
Exception' to sympy/{submodule}/__init__.py and observe the traceback from
running 'from sympy import *'.""")

    print("No unexpected submodules were imported with 'from sympy import *'")

if __name__ == '__main__':
    test_submodule_imports()
