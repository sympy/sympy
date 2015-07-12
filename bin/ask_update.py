#!/usr/bin/env python

""" Update the ask_generated.py file

This must be run each time known_facts is changed

Should be run from sympy root directory

$ python bin/ask_update.py
"""

from __future__ import division, print_function

# hook in-tree SymPy into Python path, if possible
import os
import sys

isympy_path = os.path.abspath(__file__)
isympy_dir = os.path.dirname(isympy_path)
sympy_top = os.path.split(isympy_dir)[0]
sympy_dir = os.path.join(sympy_top, 'sympy')

if os.path.isdir(sympy_dir):
    sys.path.insert(0, sympy_top)

from sympy.assumptions.ask import (compute_known_facts,
        get_known_facts, get_known_facts_keys)

with open('sympy/assumptions/ask_generated.py', 'w') as f:
    code = compute_known_facts(get_known_facts(), get_known_facts_keys())
    f.write(code)
