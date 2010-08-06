"""Module for querying SymPy objects about assumptions."""
import inspect
import copy
from sympy.core import sympify
from sympy.utilities.source import get_class
from sympy.assumptions import global_assumptions, Assume, Predicate
from sympy.assumptions.assume import eliminate_assume
from sympy.logic.boolalg import to_cnf, conjuncts, disjuncts, \
    And, Not, Or, Implies, Equivalent, to_int_repr
from sympy.logic.inference import literal_symbol, satisfiable

class Q:
    """Supported ask keys."""
    bounded = Predicate('bounded')
    commutative = Predicate('commutative')
    complex = Predicate('complex')
    composite = Predicate('composite')
    even = Predicate('even')
    extended_real = Predicate('extended_real')
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


def eval_predicate(predicate, expr, assumptions=True):
    """
    Evaluate predicate(expr) under the given assumptions.

    This uses only direct resolution methods, not logical inference.
    """
    res, _res = None, None
    mro = inspect.getmro(type(expr))
    for handler in predicate.handlers:
        cls = get_class(handler)
        for subclass in mro:
            try:
                eval = getattr(cls, subclass.__name__)
            except AttributeError:
                continue
            res = eval(expr, assumptions)
            if _res is None:
                _res = res
            elif res is None:
                # since first resolutor was conclusive, we keep that value
                res = _res
            else:
                # only check consistency if both resolutors have concluded
                if _res != res:
                    raise ValueError('incompatible resolutors')
            break
    return res


def ask(expr, key, assumptions=True, context=global_assumptions):
    """
    Method for inferring properties about objects.

    **Syntax**

        * ask(expression, key)

        * ask(expression, key, assumptions)

            where expression is any SymPy expression

    **Examples**
        >>> from sympy import ask, Q, Assume, pi
        >>> from sympy.abc import x, y
        >>> ask(pi, Q.rational)
        False
        >>> ask(x*y, Q.even, Assume(x, Q.even) & Assume(y, Q.integer))
        True
        >>> ask(x*y, Q.prime, Assume(x, Q.integer) &  Assume(y, Q.integer))
        False

    **Remarks**
        Relations in assumptions are not implemented (yet), so the following
        will not give a meaningful result.
        >> ask(x, positive=True, Assume(x>0))
        It is however a work in progress and should be available before
        the official release

    """
    expr = sympify(expr)
    if type(key) is not Predicate:
        key = getattr(Q, str(key))
    assumptions = And(assumptions, And(*context))

    # direct resolution method, no logic
    res = eval_predicate(key, expr, assumptions)
    if res is not None:
        return res

    # use logic inference
    if assumptions is True:
        return

    if not expr.is_Atom:
        return

    assumptions = eliminate_assume(assumptions, expr)
    if assumptions is None:
        return

    # If it's not consistent with the assumptions, then it can't be true
    if not satisfiable(And(known_facts_cnf, assumptions, key)):
        return False

    # If the negation is unsatisfiable, it is entailed
    if not satisfiable(And(known_facts_cnf, assumptions, Not(key))):
        return True

    # Otherwise, we don't have enough information to conclude one way or the other
    return None

def register_handler(key, handler):
    """Register a handler in the ask system. key must be a string and handler a
    class inheriting from AskHandler.

        >>> from sympy.assumptions import register_handler, ask
        >>> from sympy.assumptions.handlers import AskHandler
        >>> class MersenneHandler(AskHandler):
        ...     # Mersenne numbers are in the form 2**n + 1, n integer
        ...     @staticmethod
        ...     def Integer(expr, assumptions):
        ...         import math
        ...         return ask(math.log(expr + 1, 2), 'integer')
        >>> register_handler('mersenne', MersenneHandler)
        >>> ask(7, 'mersenne')
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
    # Compute the known facts in CNF form for logical inference
    fact_string = " -{ Known facts in CNF }-\n"
    cnf = to_cnf(known_facts)
    fact_string += "known_facts_cnf = And( \\\n   ",
    fact_string += ", \\\n    ".join(map(str, cnf.args))
    fact_string += "\n)\n"

    # Compute the quick lookup for single facts
    from sympy.abc import x
    mapping = {}
    for key in known_facts_keys:
        mapping[key] = set([])
        for other_key in known_facts_keys:
            if other_key != key:
                if ask(x, other_key, Assume(x, key)):
                    mapping[key].add(other_key)
    fact_string += "\n\n -{ Known facts in compressed sets }-\n"
    fact_string += "known_facts_dict = { \\\n   ",
    fact_string += ", \\\n    ".join(["%s: %s" % item for item in mapping.items()])
    fact_string += "\n}\n"
    return fact_string

# handlers_dict tells us what ask handler we should use
# for a particular key
_handlers_dict = {
    'bounded'        : ['sympy.assumptions.handlers.calculus.AskBoundedHandler'],
    'commutative'    : ['sympy.assumptions.handlers.AskCommutativeHandler'],
    'complex'        : ['sympy.assumptions.handlers.sets.AskComplexHandler'],
    'composite'      : ['sympy.assumptions.handlers.ntheory.AskCompositeHandler'],
    'even'           : ['sympy.assumptions.handlers.ntheory.AskEvenHandler'],
    'extended_real'  : ['sympy.assumptions.handlers.sets.AskExtendedRealHandler'],
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
    Equivalent(Q.even, Q.integer & ~Q.odd),
    Equivalent(Q.extended_real, Q.real | Q.infinity),
    Equivalent(Q.odd, Q.integer & ~Q.even),
    Equivalent(Q.prime, Q.integer & Q.positive & ~Q.composite),
    Implies   (Q.integer, Q.rational),
    Implies   (Q.imaginary, Q.complex & ~Q.real),
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

known_facts_cnf = And( \
    Or(Not(Q.integer), Q.even, Q.odd), \
    Or(Not(Q.extended_real), Q.infinity, Q.real), \
    Or(Not(Q.real), Q.irrational, Q.rational), \
    Or(Not(Q.real), Q.complex), \
    Or(Not(Q.integer), Not(Q.positive), Q.composite, Q.prime), \
    Or(Not(Q.integer), Q.rational), \
    Or(Not(Q.imaginary), Q.complex), \
    Or(Not(Q.even), Q.integer), \
    Or(Not(Q.positive), Q.nonzero), \
    Or(Not(Q.nonzero), Q.negative, Q.positive), \
    Or(Not(Q.prime), Q.positive), \
    Or(Not(Q.rational), Q.real), \
    Or(Not(Q.imaginary), Not(Q.real)), \
    Or(Not(Q.odd), Q.integer), \
    Or(Not(Q.real), Q.extended_real), \
    Or(Not(Q.composite), Not(Q.prime)), \
    Or(Not(Q.negative), Q.nonzero), \
    Or(Not(Q.negative), Not(Q.positive)), \
    Or(Not(Q.prime), Q.integer), \
    Or(Not(Q.even), Not(Q.odd)), \
    Or(Not(Q.nonzero), Q.real), \
    Or(Not(Q.irrational), Q.real), \
    Or(Not(Q.irrational), Not(Q.rational)), \
    Or(Not(Q.infinity), Q.extended_real)
)

known_facts_dict = { \
    Q.is_true: set([Q.commutative]), \
    Q.complex: set([Q.commutative]), \
    Q.odd: set([Q.complex, Q.real, Q.rational, Q.extended_real, Q.integer, Q.commutative]), \
    Q.positive: set([Q.real, Q.commutative, Q.complex, Q.extended_real, Q.nonzero]), \
    Q.real: set([Q.extended_real, Q.commutative, Q.complex]), \
    Q.composite: set([Q.commutative]), \
    Q.bounded: set([Q.commutative]), \
    Q.prime: set([Q.complex, Q.positive, Q.nonzero, Q.real, Q.rational, Q.extended_real, Q.integer, Q.commutative]), \
    Q.infinitesimal: set([Q.commutative]), \
    Q.even: set([Q.complex, Q.real, Q.rational, Q.extended_real, Q.integer, Q.commutative]), \
    Q.negative: set([Q.real, Q.commutative, Q.complex, Q.extended_real, Q.nonzero]), \
    Q.rational: set([Q.real, Q.commutative, Q.complex, Q.extended_real]), \
    Q.extended_real: set([Q.commutative]), \
    Q.nonzero: set([Q.real, Q.commutative, Q.complex, Q.extended_real]), \
    Q.integer: set([Q.real, Q.commutative, Q.rational, Q.complex, Q.extended_real]), \
    Q.irrational: set([Q.real, Q.commutative, Q.complex, Q.extended_real]), \
    Q.commutative: set([]), \
    Q.infinity: set([Q.commutative, Q.extended_real]), \
    Q.algebraic: set([Q.commutative]), \
    Q.imaginary: set([Q.commutative, Q.complex])
}
