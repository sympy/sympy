import inspect
import copy
from sympy.core import Symbol, sympify
from sympy.utilities.source import get_class
from sympy.assumptions import global_assumptions
from sympy.assumptions.assume import eliminate_assume
from sympy.logic.boolalg import to_cnf, conjuncts, disjuncts, \
    Equivalent, And, Or, Not
from sympy.logic.inference import literal_symbol
from sympy.logic.algorithms.dpll import dpll_int_repr

class Q:
    """Supported ask keys"""
    bounded = 'bounded'
    commutative = 'commutative'
    complex = 'complex'
    composite = 'composite'
    even = 'even'
    extended_real = 'extended_real'
    imaginary = 'imaginary'
    infinitesimal = 'infinitesimal'
    infinity = 'infinity'
    integer = 'integer'
    irrational = 'irrational'
    rational = 'rational'
    negative = 'negative'
    nonzero = 'nonzero'
    positive = 'positive'
    prime = 'prime'
    real = 'real'
    odd = 'odd'

# TODO: maybe this should be moved to another file?
def ask(expr, key, assumptions=True):
    """
    Method for inferring properties about objects.

    Syntax

        * ask(expression, key)

        * ask(expression, key, assumptions)

            where expression is any SymPy expression

    Examples
        >>> from sympy import ask, Q, Assume, pi
        >>> from sympy.abc import x, y
        >>> ask(pi, Q.rational)
        False
        >>> ask(x*y, Q.even, Assume(x, Q.even) & Assume(y, Q.integer))
        True
        >>> ask(x*y, Q.prime, Assume(x, Q.integer) &  Assume(y, Q.integer))
        False

    Remarks
        Relations in assumptions are not implemented (yet), so the following
        will not give a meaningful result.
        >> ask(x, positive=True, Assume(x>0))
        It is however a work in progress and should be available before
        the official release

    """
    expr = sympify(expr)
    assumptions = And(assumptions, And(*global_assumptions))

    # direct resolution method, no logic
    resolutors = []
    for handler in handlers_dict[key]:
        resolutors.append( get_class(handler) )
    res, _res = None, None
    mro = inspect.getmro(type(expr))
    for handler in resolutors:
        for subclass in mro:
            if hasattr(handler, subclass.__name__):
                res = getattr(handler, subclass.__name__)(expr, assumptions)
                if _res is None: _res = res
                elif _res != res: raise ValueError, 'incompatible resolutors'
                break
    if res is not None:
        return res

    if assumptions is True: return

    # use logic inference
    if not expr.is_Atom: return
    clauses = copy.deepcopy(known_facts_compiled)

    assumptions = conjuncts(to_cnf(assumptions))
    # add assumptions to the knowledge base
    for assump in assumptions:
        conj = eliminate_assume(assump, symbol=expr)
        if conj:
            out = []
            for sym in conjuncts(to_cnf(conj)):
                lit, pos = literal_symbol(sym), type(sym) is not Not
                if pos:
                    out.extend([known_facts_keys.index(str(l))+1 for l in disjuncts(lit)])
                else:
                    out.extend([-(known_facts_keys.index(str(l))+1) for l in disjuncts(lit)])
            clauses.append(out)

    n = len(known_facts_keys)
    clauses.append([known_facts_keys.index(key)+1])
    if not dpll_int_repr(clauses, range(1, n+1), {}):
        return False
    clauses[-1][0] = -clauses[-1][0]
    if not dpll_int_repr(clauses, range(1, n+1), {}):
        # if the negation is satisfiable, it is entailed
        return True
    del clauses


def register_handler(key, handler):
    """Register a handler in the ask system. key must be a string and handler a
    class inheriting from AskHandler.

        >>> from sympy.queries import register_handler, ask
        >>> from sympy.queries.handlers import AskHandler
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
    if key in handlers_dict:
        handlers_dict[key].append(handler)
    else:
        handlers_dict.update({key: [handler]})

def remove_handler(key, handler):
    """Removes a handler from the ask system. Same syntax as register_handler"""
    handlers_dict[key].remove(handler)

# handlers_dict tells us what ask handler we should use
# for a particular key
handlers_dict = {
    'bounded'        : ['sympy.queries.handlers.calculus.AskBoundedHandler'],
    'commutative'    : ['sympy.queries.handlers.AskCommutativeHandler'],
    'complex'        : ['sympy.queries.handlers.sets.AskComplexHandler'],
    'composite'      : ['sympy.queries.handlers.ntheory.AskCompositeHandler'],
    'even'           : ['sympy.queries.handlers.ntheory.AskEvenHandler'],
    'extended_real'  : ['sympy.queries.handlers.sets.AskExtendedRealHandler'],
    'imaginary'      : ['sympy.queries.handlers.sets.AskImaginaryHandler'],
    'infinitesimal'  : ['sympy.queries.handlers.calculus.AskInfinitesimalHandler'],
    'integer'        : ['sympy.queries.handlers.sets.AskIntegerHandler'],
    'irrational'     : ['sympy.queries.handlers.sets.AskIrrationalHandler'],
    'rational'       : ['sympy.queries.handlers.sets.AskRationalHandler'],
    'negative'       : ['sympy.queries.handlers.order.AskNegativeHandler'],
    'nonzero'        : ['sympy.queries.handlers.order.AskNonZeroHandler'],
    'positive'       : ['sympy.queries.handlers.order.AskPositiveHandler'],
    'prime'          : ['sympy.queries.handlers.ntheory.AskPrimeHandler'],
    'real'           : ['sympy.queries.handlers.sets.AskRealHandler'],
    'odd'            : ['sympy.queries.handlers.ntheory.AskOddHandler'],
}

known_facts_keys = []

for k in Q.__dict__.keys():
    if k.startswith('__'): continue
    known_facts_keys.append(k)

"""
known_facts_compiled gets generated from the above. known_facts_compiled
is in CNF form and efficient integer representation.

known_facts = [
    Implies   (real, complex),
    Equivalent(even, integer & ~odd),
    Equivalent(extended_real, real | infinity),
    Equivalent(odd, integer & ~even),
    Equivalent(prime, integer & positive & ~composite),
    Implies   (integer, rational),
    Implies   (imaginary, complex & ~real),
    Equivalent(negative, nonzero & ~positive),
    Equivalent(positive, nonzero & ~negative),
    Equivalent(rational, real & ~irrational),
    Equivalent(real, rational | irrational),
    Implies   (nonzero, real),
    Equivalent(nonzero, positive | negative),
]

To generate known_facts_compiled, use the following script.

from sympy import var
from sympy.queries import Q
from sympy.logic.boolalg import *

syms = []

for k in Q.__dict__.keys():
    if k.startswith('__'): continue
    syms.append(var(k))
print syms
out = []

for _c in known_facts:
    _c = conjuncts(to_cnf(_c))
    c = to_int_repr(_c, syms)
    for l in c:
        out.append(l)
print out
"""
known_facts_compiled = [[11, -14], [15, -1], [-1, -17], [1, 17, -15], [3, -4], [3, -14], [4, 14, -3], [15, -17], [-1, -17], [1, 17, -15], [15, -10], [7, -10], [-6, -10], [6, 10, -15, -7], [13, -15], [11, -16], [-16, -14], [2, -9], [-9, -7], [9, 7, -2], [2, -7], [-9, -7], [9, 7, -2], [14, -13], [-18, -13], [18, 13, -14], [14, -18], [14, -13], [18, 13, -14], [14, -2], [2, -9], [2, -7], [9, 7, -2]]





