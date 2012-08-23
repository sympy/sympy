"""Module for querying SymPy objects about assumptions."""
from sympy.core import sympify
from sympy.logic.boolalg import to_cnf, And, Not, Or, Implies, Equivalent
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
    negative = Predicate('negative')
    nonzero = Predicate('nonzero')
    positive = Predicate('positive')
    prime = Predicate('prime')
    real = Predicate('real')
    odd = Predicate('odd')
    is_true = Predicate('is_true')



def _extract_facts(expr, symbol):
    """
    Helper for ask().

    Extracts the facts relevant to the symbol from an assumption.
    Returns None if there is nothing to extract.
    """
    if not expr.has(symbol):
        return None
    if isinstance(expr, AppliedPredicate):
        return expr.func
    return expr.func(*filter(lambda x: x is not None,
                [_extract_facts(arg, symbol) for arg in expr.args]))

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
    assumptions = And(assumptions, And(*context))
    if isinstance(proposition, AppliedPredicate):
        key, expr = proposition.func, sympify(proposition.arg)
    else:
        key, expr = Q.is_true, sympify(proposition)

    # direct resolution method, no logic
    res = key(expr)._eval_ask(assumptions)
    if res is not None:
        return res

    if assumptions is True:
        return

    if not expr.is_Atom:
        return

    local_facts = _extract_facts(assumptions, expr)
    if local_facts is None or local_facts is True:
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
    return ask_full_inference(key, local_facts)


def ask_full_inference(proposition, assumptions):
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

def compute_known_facts():
    """Compute the various forms of knowledge compilation used by the
    assumptions system.
    """
    # Compute the known facts in CNF form for logical inference
    fact_string = "# -{ Known facts in CNF }-\n"
    cnf = to_cnf(known_facts)
    fact_string += "known_facts_cnf = And(\n    "
    fact_string += ",\n    ".join(map(str, cnf.args))
    fact_string += "\n)\n"

    # Compute the quick lookup for single facts
    mapping = {}
    for key in known_facts_keys:
        mapping[key] = set([key])
        for other_key in known_facts_keys:
            if other_key != key:
                if ask_full_inference(other_key, key):
                    mapping[key].add(other_key)
    fact_string += "\n# -{ Known facts in compressed sets }-\n"
    fact_string += "known_facts_dict = {\n    "
    fact_string += ",\n    ".join(["%s: %s" % item for item in mapping.items()])
    fact_string += "\n}\n"
    return fact_string

# handlers_dict tells us what ask handler we should use
# for a particular key
_handlers_dict = {
    'antihermitian'  : ['sympy.assumptions.handlers.sets.AskAntiHermitianHandler'],
    'bounded'        : ['sympy.assumptions.handlers.calculus.AskBoundedHandler'],
    'commutative'    : ['sympy.assumptions.handlers.AskCommutativeHandler'],
    'complex'        : ['sympy.assumptions.handlers.sets.AskComplexHandler'],
    'composite'      : ['sympy.assumptions.handlers.ntheory.AskCompositeHandler'],
    'even'           : ['sympy.assumptions.handlers.ntheory.AskEvenHandler'],
    'extended_real'  : ['sympy.assumptions.handlers.sets.AskExtendedRealHandler'],
    'hermitian'      : ['sympy.assumptions.handlers.sets.AskHermitianHandler'],
    'imaginary'      : ['sympy.assumptions.handlers.sets.AskImaginaryHandler'],
    'infinitesimal'  : ['sympy.assumptions.handlers.calculus.AskInfinitesimalHandler'],
    'integer'        : ['sympy.assumptions.handlers.sets.AskIntegerHandler'],
    'irrational'     : ['sympy.assumptions.handlers.sets.AskIrrationalHandler'],
    'rational'       : ['sympy.assumptions.handlers.sets.AskRationalHandler'],
    'negative'       : ['sympy.assumptions.handlers.order.AskNegativeHandler'],
    'nonzero'        : ['sympy.assumptions.handlers.order.AskNonZeroHandler'],
    'positive'       : ['sympy.assumptions.handlers.order.AskPositiveHandler'],
    'prime'          : ['sympy.assumptions.handlers.ntheory.AskPrimeHandler'],
    'real'           : ['sympy.assumptions.handlers.sets.AskRealHandler'],
    'odd'            : ['sympy.assumptions.handlers.ntheory.AskOddHandler'],
    'algebraic'      : ['sympy.assumptions.handlers.sets.AskAlgebraicHandler'],
    'is_true'        : ['sympy.assumptions.handlers.TautologicalHandler']
}
for name, value in _handlers_dict.iteritems():
    register_handler(name, value[0])


known_facts_keys = [getattr(Q, attr) for attr in Q.__dict__ \
                                                if not attr.startswith('__')]
known_facts = And(
    Implies   (Q.real, Q.complex),
    Implies   (Q.real, Q.hermitian),
    Equivalent(Q.even, Q.integer & ~Q.odd),
    Equivalent(Q.extended_real, Q.real | Q.infinity),
    Equivalent(Q.odd, Q.integer & ~Q.even),
    Equivalent(Q.prime, Q.integer & Q.positive & ~Q.composite),
    Implies   (Q.integer, Q.rational),
    Implies   (Q.imaginary, Q.complex & ~Q.real),
    Implies   (Q.imaginary, Q.antihermitian),
    Implies   (Q.antihermitian, ~Q.hermitian),
    Equivalent(Q.negative, Q.nonzero & ~Q.positive),
    Equivalent(Q.positive, Q.nonzero & ~Q.negative),
    Equivalent(Q.rational, Q.real & ~Q.irrational),
    Equivalent(Q.real, Q.rational | Q.irrational),
    Implies   (Q.nonzero, Q.real),
    Equivalent(Q.nonzero, Q.positive | Q.negative)
)

################################################################################
# Note: The following facts are generated by the compute_known_facts function. #
################################################################################
# -{ Known facts in CNF }-
known_facts_cnf = And(
    Or(Not(Q.integer), Q.even, Q.odd),
    Or(Not(Q.extended_real), Q.real, Q.infinity),
    Or(Not(Q.real), Q.irrational, Q.rational),
    Or(Not(Q.real), Q.complex),
    Or(Not(Q.integer), Not(Q.positive), Q.prime, Q.composite),
    Or(Not(Q.imaginary), Q.antihermitian),
    Or(Not(Q.integer), Q.rational),
    Or(Not(Q.real), Q.hermitian),
    Or(Not(Q.imaginary), Q.complex),
    Or(Not(Q.even), Q.integer),
    Or(Not(Q.positive), Q.nonzero),
    Or(Not(Q.nonzero), Q.negative, Q.positive),
    Or(Not(Q.prime), Q.positive),
    Or(Not(Q.rational), Q.real),
    Or(Not(Q.real), Not(Q.imaginary)),
    Or(Not(Q.odd), Q.integer),
    Or(Not(Q.real), Q.extended_real),
    Or(Not(Q.composite), Not(Q.prime)),
    Or(Not(Q.negative), Q.nonzero),
    Or(Not(Q.positive), Not(Q.negative)),
    Or(Not(Q.prime), Q.integer),
    Or(Not(Q.even), Not(Q.odd)),
    Or(Not(Q.nonzero), Q.real),
    Or(Not(Q.irrational), Q.real),
    Or(Not(Q.rational), Not(Q.irrational)),
    Or(Not(Q.infinity), Q.extended_real),
    Or(Not(Q.antihermitian), Not(Q.hermitian))
)

# -{ Known facts in compressed sets }-
known_facts_dict = {
    Q.odd: set([Q.complex, Q.odd, Q.hermitian, Q.real, Q.rational, Q.extended_real, Q.integer]),
    Q.antihermitian: set([Q.antihermitian]),
    Q.infinitesimal: set([Q.infinitesimal]),
    Q.hermitian: set([Q.hermitian]),
    Q.bounded: set([Q.bounded]),
    Q.even: set([Q.complex, Q.real, Q.hermitian, Q.even, Q.rational, Q.extended_real, Q.integer]),
    Q.algebraic: set([Q.algebraic]),
    Q.is_true: set([Q.is_true]),
    Q.real: set([Q.real, Q.complex, Q.extended_real, Q.hermitian]),
    Q.rational: set([Q.real, Q.rational, Q.complex, Q.extended_real, Q.hermitian]),
    Q.extended_real: set([Q.extended_real]),
    Q.integer: set([Q.complex, Q.hermitian, Q.real, Q.rational, Q.extended_real, Q.integer]),
    Q.commutative: set([Q.commutative]),
    Q.infinity: set([Q.extended_real, Q.infinity]),
    Q.complex: set([Q.complex]),
    Q.positive: set([Q.complex, Q.positive, Q.nonzero, Q.hermitian, Q.real, Q.extended_real]),
    Q.composite: set([Q.composite]),
    Q.prime: set([Q.complex, Q.positive, Q.real, Q.hermitian, Q.prime, Q.rational, Q.extended_real, Q.nonzero, Q.integer]),
    Q.negative: set([Q.complex, Q.nonzero, Q.hermitian, Q.real, Q.negative, Q.extended_real]),
    Q.nonzero: set([Q.nonzero, Q.complex, Q.extended_real, Q.real, Q.hermitian]),
    Q.irrational: set([Q.real, Q.irrational, Q.complex, Q.extended_real, Q.hermitian]),
    Q.imaginary: set([Q.antihermitian, Q.complex, Q.imaginary])
}
