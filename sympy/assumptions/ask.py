"""Module for querying SymPy objects about assumptions."""
from __future__ import print_function, division

from sympy.core import sympify
from sympy.logic.boolalg import (to_cnf, And, Not, Or, Implies, Equivalent,
    BooleanFunction, BooleanAtom)
from sympy.logic.inference import satisfiable
from sympy.assumptions.assume import (global_assumptions, Predicate,
        AppliedPredicate)


class Q:
    """Supported ask keys."""
    antihermitian = Predicate('antihermitian')
    bounded = Predicate('bounded')
    commutative = Predicate('commutative')
    complex = Predicate('complex')
    composite = Predicate('composite')
    even = Predicate('even')
    extended_real = Predicate('extended_real')
    hermitian = Predicate('hermitian')
    imaginary = Predicate('imaginary')
    infinitesimal = Predicate('infinitesimal')
    infinity = Predicate('infinity')
    integer = Predicate('integer')
    irrational = Predicate('irrational')
    rational = Predicate('rational')
    algebraic = Predicate('algebraic')
    transcendental = Predicate('transcendental')
    negative = Predicate('negative')
    nonzero = Predicate('nonzero')
    positive = Predicate('positive')
    prime = Predicate('prime')
    real = Predicate('real')
    odd = Predicate('odd')
    is_true = Predicate('is_true')
    nonpositive = Predicate('nonpositive')
    nonnegative = Predicate('nonnegative')
    zero = Predicate('zero')

    symmetric = Predicate('symmetric')
    invertible = Predicate('invertible')
    singular = Predicate('singular')
    orthogonal = Predicate('orthogonal')
    unitary = Predicate('unitary')
    normal = Predicate('normal')
    positive_definite = Predicate('positive_definite')
    upper_triangular = Predicate('upper_triangular')
    lower_triangular = Predicate('lower_triangular')
    diagonal = Predicate('diagonal')
    triangular = Predicate('triangular')
    unit_triangular = Predicate('unit_triangular')
    fullrank = Predicate('fullrank')
    square = Predicate('square')
    real_elements = Predicate('real_elements')
    complex_elements = Predicate('complex_elements')
    integer_elements = Predicate('integer_elements')


def _extract_facts(expr, symbol):
    """
    Helper for ask().

    Extracts the facts relevant to the symbol from an assumption.
    Returns None if there is nothing to extract.
    """
    if isinstance(expr, bool):
        return
    if not expr.has(symbol):
        return
    if isinstance(expr, AppliedPredicate):
        if expr.arg == symbol:
            return expr.func
        else:
            return
    if isinstance(expr, Not) and expr.args[0].func in (And, Or):
        cls = Or if expr.args[0] == And else And
        expr = cls(*[~arg for arg in expr.args[0].args])
    args = [_extract_facts(arg, symbol) for arg in expr.args]
    if isinstance(expr, And):
        args = [x for x in args if x is not None]
        if args:
            return expr.func(*args)
    if args and all(x is not None for x in args):
        return expr.func(*args)


def ask(proposition, assumptions=True, context=global_assumptions):
    """
    Method for inferring properties about objects.

    **Syntax**

        * ask(proposition)

        * ask(proposition, assumptions)

            where ``proposition`` is any boolean expression

    Examples
    ========

    >>> from sympy import ask, Q, pi
    >>> from sympy.abc import x, y
    >>> ask(Q.rational(pi))
    False
    >>> ask(Q.even(x*y), Q.even(x) & Q.integer(y))
    True
    >>> ask(Q.prime(x*y), Q.integer(x) &  Q.integer(y))
    False

    **Remarks**
        Relations in assumptions are not implemented (yet), so the following
        will not give a meaningful result.

        >>> ask(Q.positive(x), Q.is_true(x > 0)) # doctest: +SKIP

        It is however a work in progress.

    """
    if not isinstance(proposition, (BooleanFunction, AppliedPredicate, bool, BooleanAtom)):
        raise TypeError("proposition must be a valid logical expression")

    if not isinstance(assumptions, (BooleanFunction, AppliedPredicate, bool, BooleanAtom)):
        raise TypeError("assumptions must be a valid logical expression")

    if isinstance(proposition, AppliedPredicate):
        key, expr = proposition.func, sympify(proposition.arg)
    else:
        key, expr = Q.is_true, sympify(proposition)

    assumptions = And(assumptions, And(*context))
    assumptions = to_cnf(assumptions)

    local_facts = _extract_facts(assumptions, expr)

    if local_facts and satisfiable(And(local_facts, known_facts_cnf)) is False:
        raise ValueError("inconsistent assumptions %s" % assumptions)

    # direct resolution method, no logic
    res = key(expr)._eval_ask(assumptions)
    if res is not None:
        return bool(res)

    if assumptions is True:
        return

    if local_facts is None:
        return

    # See if there's a straight-forward conclusion we can make for the inference
    if local_facts.is_Atom:
        if key in known_facts_dict[local_facts]:
            return True
        if Not(key) in known_facts_dict[local_facts]:
            return False
    elif local_facts.func is And and all(k in known_facts_dict for k in local_facts.args):
        for assum in local_facts.args:
            if assum.is_Atom:
                if key in known_facts_dict[assum]:
                    return True
                if Not(key) in known_facts_dict[assum]:
                    return False
            elif assum.func is Not and assum.args[0].is_Atom:
                if key in known_facts_dict[assum]:
                    return False
                if Not(key) in known_facts_dict[assum]:
                    return True
    elif (isinstance(key, Predicate) and
            local_facts.func is Not and local_facts.args[0].is_Atom):
        if local_facts.args[0] in known_facts_dict[key]:
            return False

    # Failing all else, we do a full logical inference
    return ask_full_inference(key, local_facts, known_facts_cnf)


def ask_full_inference(proposition, assumptions, known_facts_cnf):
    """
    Method for inferring properties about objects.

    """
    if not satisfiable(And(known_facts_cnf, assumptions, proposition)):
        return False
    if not satisfiable(And(known_facts_cnf, assumptions, Not(proposition))):
        return True
    return None


def register_handler(key, handler):
    """
    Register a handler in the ask system. key must be a string and handler a
    class inheriting from AskHandler::

        >>> from sympy.assumptions import register_handler, ask, Q
        >>> from sympy.assumptions.handlers import AskHandler
        >>> class MersenneHandler(AskHandler):
        ...     # Mersenne numbers are in the form 2**n + 1, n integer
        ...     @staticmethod
        ...     def Integer(expr, assumptions):
        ...         import math
        ...         return ask(Q.integer(math.log(expr + 1, 2)))
        >>> register_handler('mersenne', MersenneHandler)
        >>> ask(Q.mersenne(7))
        True

    """
    if type(key) is Predicate:
        key = key.name
    try:
        getattr(Q, key).add_handler(handler)
    except AttributeError:
        setattr(Q, key, Predicate(key, handlers=[handler]))


def remove_handler(key, handler):
    """Removes a handler from the ask system. Same syntax as register_handler"""
    if type(key) is Predicate:
        key = key.name
    getattr(Q, key).remove_handler(handler)


def single_fact_lookup(known_facts_keys, known_facts_cnf):
    # Compute the quick lookup for single facts
    mapping = {}
    for key in known_facts_keys:
        mapping[key] = set([key])
        for other_key in known_facts_keys:
            if other_key != key:
                if ask_full_inference(other_key, key, known_facts_cnf):
                    mapping[key].add(other_key)
    return mapping


def compute_known_facts(known_facts, known_facts_keys):
    """Compute the various forms of knowledge compilation used by the
    assumptions system.

    This function is typically applied to the variables
    ``known_facts`` and ``known_facts_keys`` defined at the bottom of
    this file.
    """
    from textwrap import dedent, wrap

    fact_string = dedent('''\
    """
    The contents of this file are the return value of
    ``sympy.assumptions.ask.compute_known_facts``.  Do NOT manually
    edit this file.  Instead, run ./bin/ask_update.py.
    """

    from sympy.logic.boolalg import And, Not, Or
    from sympy.assumptions.ask import Q

    # -{ Known facts in CNF }-
    known_facts_cnf = And(
        %s
    )

    # -{ Known facts in compressed sets }-
    known_facts_dict = {
        %s
    }
    ''')
    # Compute the known facts in CNF form for logical inference
    LINE = ",\n    "
    HANG = ' '*8
    cnf = to_cnf(known_facts)
    c = LINE.join([str(a) for a in cnf.args])
    mapping = single_fact_lookup(known_facts_keys, cnf)
    items = sorted(mapping.items(), key=str)
    keys = [str(i[0]) for i in items]
    values = ['set(%s)' % sorted(i[1], key=str) for i in items]
    m = LINE.join(['\n'.join(
        wrap("%s: %s" % (k, v),
            subsequent_indent=HANG,
            break_long_words=False))
        for k, v in zip(keys, values)]) + ','
    return fact_string % (c, m)

# handlers tells us what ask handler we should use
# for a particular key
_val_template = 'sympy.assumptions.handlers.%s'
_handlers = [
    ("antihermitian",     "sets.AskAntiHermitianHandler"),
    ("bounded",           "calculus.AskBoundedHandler"),
    ("commutative",       "AskCommutativeHandler"),
    ("complex",           "sets.AskComplexHandler"),
    ("composite",         "ntheory.AskCompositeHandler"),
    ("even",              "ntheory.AskEvenHandler"),
    ("extended_real",     "sets.AskExtendedRealHandler"),
    ("hermitian",         "sets.AskHermitianHandler"),
    ("imaginary",         "sets.AskImaginaryHandler"),
    ("infinitesimal",     "calculus.AskInfinitesimalHandler"),
    ("integer",           "sets.AskIntegerHandler"),
    ("irrational",        "sets.AskIrrationalHandler"),
    ("rational",          "sets.AskRationalHandler"),
    ("negative",          "order.AskNegativeHandler"),
    ("nonzero",           "order.AskNonZeroHandler"),
    ("nonpositive",       "order.AskNonPositiveHandler"),
    ("nonnegative",       "order.AskNonNegativeHandler"),
    ("zero",              "order.AskZeroHandler"),
    ("positive",          "order.AskPositiveHandler"),
    ("prime",             "ntheory.AskPrimeHandler"),
    ("real",              "sets.AskRealHandler"),
    ("odd",               "ntheory.AskOddHandler"),
    ("algebraic",         "sets.AskAlgebraicHandler"),
    ("is_true",           "common.TautologicalHandler"),
    ("symmetric",         "matrices.AskSymmetricHandler"),
    ("invertible",        "matrices.AskInvertibleHandler"),
    ("orthogonal",        "matrices.AskOrthogonalHandler"),
    ("unitary",           "matrices.AskUnitaryHandler"),
    ("positive_definite", "matrices.AskPositiveDefiniteHandler"),
    ("upper_triangular",  "matrices.AskUpperTriangularHandler"),
    ("lower_triangular",  "matrices.AskLowerTriangularHandler"),
    ("diagonal",          "matrices.AskDiagonalHandler"),
    ("fullrank",          "matrices.AskFullRankHandler"),
    ("square",            "matrices.AskSquareHandler"),
    ("integer_elements",  "matrices.AskIntegerElementsHandler"),
    ("real_elements",     "matrices.AskRealElementsHandler"),
    ("complex_elements",  "matrices.AskComplexElementsHandler"),
]

for name, value in _handlers:
    register_handler(name, _val_template % value)

known_facts_keys = [getattr(Q, attr) for attr in Q.__dict__
                    if not attr.startswith('__')]
known_facts = And(
    Implies(Q.real, Q.complex),
    Implies(Q.real, Q.hermitian),
    Equivalent(Q.even, Q.integer & ~Q.odd),
    Equivalent(Q.extended_real, Q.real | Q.infinity),
    Equivalent(Q.odd, Q.integer & ~Q.even),
    Equivalent(Q.prime, Q.integer & Q.positive & ~Q.composite),
    Implies(Q.integer, Q.rational),
    Implies(Q.rational, Q.algebraic),
    Implies(Q.algebraic, Q.complex),
    Equivalent(Q.transcendental, Q.complex & ~Q.algebraic),
    Implies(Q.imaginary, Q.complex & ~Q.real),
    Implies(Q.imaginary, Q.antihermitian),
    Implies(Q.antihermitian, ~Q.hermitian),
    Equivalent(Q.negative, Q.nonzero & ~Q.positive),
    Equivalent(Q.positive, Q.nonzero & ~Q.negative),
    Equivalent(Q.rational, Q.real & ~Q.irrational),
    Equivalent(Q.real, Q.rational | Q.irrational),
    Implies(Q.nonzero, Q.real),
    Equivalent(Q.nonzero, Q.positive | Q.negative),
    Equivalent(Q.nonpositive, ~Q.positive & Q.real),
    Equivalent(Q.nonnegative, ~Q.negative & Q.real),
    Equivalent(Q.zero, Q.real & ~Q.nonzero),
    Implies(Q.zero, Q.even),

    Implies(Q.orthogonal, Q.positive_definite),
    Implies(Q.orthogonal, Q.unitary),
    Implies(Q.unitary & Q.real, Q.orthogonal),
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

from sympy.assumptions.ask_generated import known_facts_dict, known_facts_cnf
