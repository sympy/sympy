#!/usr/bin/env python

""" Update the ask_generated.py file

This must be run each time known_facts is changed

Should be run from sympy root directory

$ python bin/ask_update.py
"""

from sympy.assumptions.ask import (compute_known_facts, known_facts,
        known_facts_keys)

f = open('sympy/assumptions/ask_generated.py', 'w')
code = compute_known_facts(known_facts, known_facts_keys)
f.write(code)
f.close()
