import inspect
from sympy.core import Symbol, sympify
from sympy.utilities.source import get_class
from sympy.assumptions import list_global_assumptions
from sympy.assumptions.assume import eliminate_assume
from sympy.logic.boolalg import to_cnf, conjuncts, \
    compile_rule, Equivalent, And
from sympy.logic.algorithms.dpll import dpll_satisfiable

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
    integer = 'integer'
    irrational = 'irrational'
    rational = 'rational'
    negative = 'negative'
    nonzero = 'nonzero'
    positive = 'positive'
    prime = 'prime'
    real = 'real'
    odd = 'odd'

def ask(expr, key, assumptions=[]):
    """
    Method for inferring properties about objects.

    Syntax

        * ask(expression, key)

        * ask(expression, key, assumptions)

            where expression is any SymPy expression

    Examples
        >>> from sympy import *
        >>> x, y = symbols('x y')
        >>> ask(pi, Q.rational)
        False
        >>> ask(x*y, Q.even, Assume(x, Q.even) & Assume(y, Q.integer))
        True
        >>> ask(x*y, Q.prime, Assume(x, Q.integer) &  Assume(y, Q.integer))
        False

    Remarks
        Relations in assumptions are not implemented (yet), so the following
        will not give a meaningful result.
        #>>> ask(x, positive=True, Assume(x>0))
        It is however a work in progress and should be available before
        the official release
    """
    expr = sympify(expr)

    global_assump = list_global_assumptions()
    if assumptions:
        assumptions = And(assumptions, And(*global_assump))
    elif global_assump: assumptions = And(*global_assump)
    if not isinstance(assumptions, (list, tuple)):
        assumptions = conjuncts(to_cnf(assumptions))

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

    if assumptions: pass
    else: return

    # use logic inference
    if not expr.is_Atom: return
    clauses = []
    for k, values in known_facts_dict.iteritems():
        for v in values:
            clauses.append(Equivalent(compile_rule(k), compile_rule(v)))
    result = None

    # add assumptions to the knowledge base
    for assump in assumptions:
        conj = eliminate_assume(assump, symbol=expr)
        if conj: clauses.append(conj)

    clauses.append(Symbol(key))
    # TODO: call dpll and avoid creating this object
    if not dpll_satisfiable(And(*clauses)):
        return False
    clauses[-1] = ~clauses[-1]
    if not dpll_satisfiable(And(*clauses)):
        # if the negation is satisfiable, it is entailed
        return True
    clauses.pop(-1)


def register_handler(key, handler):
    """Register a handler in the ask system. key must be a string and handler a
    class inheriting from AskHandler.

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

known_facts_dict = {
    'complex'       : ['real | complex_number_re_0'],
    'even'          : ['integer & ~odd'],
    'extended_real' : ['real | infinity'],
    'odd'           : ['integer & ~even'],
    'prime'         : ['integer & positive & ~composite'],
    'integer'       : ['rational & denominator_one'],
    'imaginary'     : ['complex & ~real'],
    'negative'      : ['real & ~positive & ~zero'],
    'nonzero'       : ['positive | negative'],
    'positive'      : ['real & ~negative & ~zero'],
    'rational'      : ['real & ~irrational'],
    'real'          : ['rational | irrational' ],
}
