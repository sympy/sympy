#!/usr/bin/env python

""" Update the ``ask_generated.py`` file.

This must be run each time ``known_facts()`` in ``assumptions.facts`` module
is changed.

This must be run each time ``_generate_assumption_rules()`` in
``sympy.core.assumptions`` module is changed.

Should be run from sympy root directory.

$ python bin/ask_update.py
"""

# hook in-tree SymPy into Python path, if possible
import os
import sys
import pprint

isympy_path = os.path.abspath(__file__)
isympy_dir = os.path.dirname(isympy_path)
sympy_top = os.path.split(isympy_dir)[0]
sympy_dir = os.path.join(sympy_top, 'sympy')

if os.path.isdir(sympy_dir):
    sys.path.insert(0, sympy_top)

from sympy.core.assumptions import _generate_assumption_rules
from sympy.assumptions.cnf import CNF, Literal
from sympy.assumptions.facts import (get_known_facts,
    generate_known_facts_dict, get_known_facts_keys, get_matrix_facts, get_number_facts)
from sympy.core import Symbol
from textwrap import dedent, wrap

def generate_code():

    LINE = ",\n        "
    HANG = ' '*8
    code_string = dedent('''\
    """
    Do NOT manually edit this file.
    Instead, run ./bin/ask_update.py.
    """

    from sympy.assumptions.ask import Q
    from sympy.assumptions.cnf import Literal
    from sympy.core.cache import cacheit

    @cacheit
    def get_all_known_facts():
        """
        Known facts between unary predicates as CNF clauses.
        """
        return {
            %s
        }

    @cacheit
    def get_all_known_matrix_facts():
        """
        Known facts between unary predicates for matrices as CNF clauses.
        """
        return {
            %s
        }

    @cacheit
    def get_all_known_number_facts():
        """
        Known facts between unary predicates for numbers as CNF clauses.
        """
        return {
            %s
        }

    @cacheit
    def get_known_facts_dict():
        """
        Logical relations between unary predicates as dictionary.

        Each key is a predicate, and item is two groups of predicates.
        First group contains the predicates which are implied by the key, and
        second group contains the predicates which are rejected by the key.

        """
        return {
            %s
        }
    ''')

    x = Symbol('x')
    fact = get_known_facts(x)
    matrix_fact = get_matrix_facts(x)
    number_fact = get_number_facts(x)

    # Generate CNF of facts between known unary predicates
    cnf = CNF.to_CNF(fact)
    all_clauses = LINE.join(sorted([
        'frozenset(('
         + ', '.join(str(Literal(lit.arg.function, lit.is_Not))
                     for lit in sorted(clause, key=str))
        + '))' for clause in cnf.clauses]))

    # Generate CNF of matrix facts
    cnf = CNF.to_CNF(matrix_fact)
    matrix_clauses = LINE.join(sorted([
        'frozenset(('
         + ', '.join(str(Literal(lit.arg.function, lit.is_Not))
                     for lit in sorted(clause, key=str))
        + '))' for clause in cnf.clauses]))

    # Generate CNF of number facts
    cnf = CNF.to_CNF(number_fact)
    number_clauses = LINE.join(sorted([
        'frozenset(('
         + ', '.join(str(Literal(lit.arg.function, lit.is_Not))
                     for lit in sorted(clause, key=str))
        + '))' for clause in cnf.clauses]))

    # Generate dictionary of facts between known unary predicates
    keys = [pred(x) for pred in get_known_facts_keys()]
    mapping = generate_known_facts_dict(keys, fact)
    items = sorted(mapping.items(), key=str)
    keys = [str(i[0]) for i in items]
    values = ['(set(%s), set(%s))' % (sorted(i[1][0], key=str),
                                      sorted(i[1][1], key=str))
              for i in items]
    m = LINE.join(['\n'.join(
        wrap("{}: {}".format(k, v),
            subsequent_indent=HANG,
            break_long_words=False))
        for k, v in zip(keys, values)]) + ','

    return code_string % (all_clauses, matrix_clauses, number_clauses, m)

with open('sympy/assumptions/ask_generated.py', 'w') as f:
    code = generate_code()
    f.write(code)


with open('sympy/core/assumptions_generated.py', 'w') as f:
    representation = _generate_assumption_rules()._to_python()

    code_string = dedent('''\
    """
    Do NOT manually edit this file.
    Instead, run ./bin/ask_update.py.
    """

    %s
    ''')

    code = code_string % (representation,)
    f.write(code)
