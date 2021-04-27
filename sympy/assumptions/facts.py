"""
Known facts in assumptions module.

This module defines the facts in ``get_known_facts()``, and supports functions
to generate the contents in ``sympy.assumptions.ask_generated`` file.
"""

from sympy.core.cache import cacheit
from sympy.assumptions import Q
from sympy.assumptions.cnf import CNF
from sympy.logic.boolalg import (to_cnf, And, Not, Implies, Equivalent)
from sympy.logic.inference import satisfiable


@cacheit
def get_composite_predicates():
    # To reduce the complexity of sat solver, these predicates never goes into facts
    # but are transformed into the combination of primitive predicates.
    return {
        Q.real : Q.negative | Q.zero | Q.positive,
        Q.integer : Q.even | Q.odd,
        Q.nonpositive : Q.negative | Q.zero,
        Q.nonzero : Q.negative | Q.positive,
        Q.nonnegative : Q.zero | Q.positive,
        Q.extended_real : Q.negative_infinite | Q.negative | Q.zero | Q.positive | Q.positive_infinite,
        Q.extended_positive: Q.positive | Q.positive_infinite,
        Q.extended_negative: Q.negative | Q.negative_infinite,
        Q.extended_nonzero: Q.negative_infinite | Q.negative | Q.positive | Q.positive_infinite,
        Q.extended_nonpositive: Q.negative_infinite | Q.negative | Q.zero,
        Q.extended_nonnegative: Q.zero | Q.positive | Q.positive_infinite,
        Q.complex : Q.algebraic | Q.transcendental
    }


@cacheit
def get_known_facts():
    # We build the facts starting with primitive predicates.
    # DO NOT include the predicates in get_composite_predicates()'s keys here!
    return And(

        # primitive predicates exclude each other
        Implies(Q.negative_infinite, ~Q.positive_infinite),
        Implies(Q.negative, ~Q.zero & ~Q.positive),
        Implies(Q.positive, ~Q.zero),

        # build real line and complex plane
        Implies(Q.negative | Q.zero | Q.positive, ~Q.imaginary),
        Implies(Q.negative | Q.zero | Q.positive | Q.imaginary, Q.algebraic | Q.transcendental),

        # other subsets of complex
        Implies(Q.transcendental, ~Q.algebraic),
        Implies(Q.irrational, ~Q.rational),
        Equivalent(Q.rational | Q.irrational, Q.negative | Q.zero | Q.positive),
        Implies(Q.rational, Q.algebraic),

        # integers
        Implies(Q.even, ~Q.odd),
        Implies(Q.even | Q.odd, Q.rational),
        Implies(Q.zero, Q.even),
        Implies(Q.composite, ~Q.prime),
        Implies(Q.composite | Q.prime, (Q.even | Q.odd) & Q.positive),
        Implies(Q.even & Q.positive & ~Q.prime, Q.composite),

        # hermitian and antihermitian
        Implies(Q.negative | Q.zero | Q.positive, Q.hermitian),
        Implies(Q.imaginary, Q.antihermitian),
        Implies(Q.zero, Q.hermitian | Q.antihermitian),

        # define finity and infinity, and build extended real line
        Implies(Q.infinite, ~Q.finite),
        Implies(Q.algebraic | Q.transcendental, Q.finite),
        Implies(Q.negative_infinite | Q.positive_infinite, Q.infinite),

        # commutativity
        Implies(Q.finite | Q.infinite, Q.commutative),

        # matrices
        Implies(Q.orthogonal, Q.positive_definite),
        Implies(Q.orthogonal, Q.unitary),
        Implies(Q.unitary & Q.real_elements, Q.orthogonal),
        Implies(Q.unitary, Q.normal),
        Implies(Q.unitary, Q.invertible),
        Implies(Q.normal, Q.square),
        Implies(Q.diagonal, Q.normal),
        Implies(Q.positive_definite, Q.invertible),
        Implies(Q.diagonal, Q.upper_triangular),
        Implies(Q.diagonal, Q.lower_triangular),
        Implies(Q.lower_triangular, Q.triangular),
        Implies(Q.upper_triangular, Q.triangular),
        Implies(Q.triangular, Q.upper_triangular | Q.lower_triangular),
        Implies(Q.upper_triangular & Q.lower_triangular, Q.diagonal),
        Implies(Q.diagonal, Q.symmetric),
        Implies(Q.unit_triangular, Q.triangular),
        Implies(Q.invertible, Q.fullrank),
        Implies(Q.invertible, Q.square),
        Implies(Q.symmetric, Q.square),
        Implies(Q.fullrank & Q.square, Q.invertible),
        Equivalent(Q.invertible, ~Q.singular),
        Implies(Q.integer_elements, Q.real_elements),
        Implies(Q.real_elements, Q.complex_elements),
    )


@cacheit
def get_known_facts_keys():
    exclude = set()
    for pred in get_composite_predicates():
        exclude.add(pred)
    for pred in [Q.eq, Q.ne, Q.gt, Q.lt, Q.ge, Q.le]:
        # sat does not support polyadic predicates yet
        exclude.add(pred)

    result = []
    for attr in Q.__class__.__dict__:
        if attr.startswith('__'):
            continue
        pred = getattr(Q, attr)
        if pred in exclude:
            continue
        result.append(pred)
    return result


def compute_known_facts(known_facts, known_facts_keys):
    """Compute the various forms of knowledge compilation used by the
    assumptions system.

    Explanation
    ===========

    This function is typically applied to the results of the ``get_known_facts``
    and ``get_known_facts_keys`` functions defined at the bottom of
    this file.
    """
    from textwrap import dedent, wrap

    fact_string = dedent('''\
    """
    The contents of this file are the return value of
    ``sympy.assumptions.ask.compute_known_facts``.

    Do NOT manually edit this file.
    Instead, run ./bin/ask_update.py.
    """

    from sympy.core.cache import cacheit
    from sympy.assumptions.cnf import Literal
    from sympy.assumptions.ask import Q

    @cacheit
    def get_all_known_facts():
        """
        Known facts as CNF clauses. Used by satask.
        """
        return {
            %s
        }

    # -{ Known facts in compressed sets }-
    @cacheit
    def get_known_facts_dict():
        """
        Logical implication as dictionary. Key implies every item in its value.
        Used for quick lookup of single facts.
        """
        return {
            %s
        }
    ''')
    # Compute the known facts in CNF form for logical inference
    LINE = ",\n        "
    HANG = ' '*8
    cnf = to_cnf(known_facts)
    cnf_ = CNF.to_CNF(known_facts)

    p = LINE.join(sorted(['frozenset((' + ', '.join(str(lit) for lit in sorted(clause, key=str)) +'))' for clause in cnf_.clauses]))
    mapping = single_fact_lookup(known_facts_keys, cnf)
    items = sorted(mapping.items(), key=str)
    keys = [str(i[0]) for i in items]
    values = ['set(%s)' % sorted(i[1], key=str) for i in items]
    m = LINE.join(['\n'.join(
        wrap("{}: {}".format(k, v),
            subsequent_indent=HANG,
            break_long_words=False))
        for k, v in zip(keys, values)]) + ','
    return fact_string % (p, m)


def single_fact_lookup(known_facts_keys, known_facts_cnf):
    # Return the dictionary for quick lookup of single fact
    mapping = {}
    for key in known_facts_keys:
        mapping[key] = {key}
        for other_key in known_facts_keys:
            if other_key != key:
                if ask_full_inference(other_key, key, known_facts_cnf):
                    mapping[key].add(other_key)
                if ask_full_inference(~other_key, key, known_facts_cnf):
                    mapping[key].add(~other_key)
    return mapping


def ask_full_inference(proposition, assumptions, known_facts_cnf):
    """
    Method for inferring properties about objects.

    """
    if not satisfiable(And(known_facts_cnf, assumptions, proposition)):
        return False
    if not satisfiable(And(known_facts_cnf, assumptions, Not(proposition))):
        return True
    return None
