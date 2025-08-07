"""
EUF SAT Task API

Module to provide a public API function to check satisfiability of formulas
in the theory of Equality with Uninterpreted Functions (EUF), using SymPy's
DPLL(T) SAT solver integrated with the EUF theory solver.

The function translates a proposition to CNF and EncodedCNF, preprocesses
predicates, and then calls the SAT solver with the EUF theory solver.

The API returns:
  - True if satisfiable,
  - False if provably unsatisfiable,
  - None if undetermined (both satisfiable and unsatisfiable assignments possible)
  
Example
-------

>>> from sympy import symbols, Eq, Function
>>> from sympy.assumptions.euf_satask import euf_satask
>>> f = Function('f')
>>> a, b, x = symbols('a b x')
>>> phi = Eq(a, b) & Eq(f(a), x)
>>> result = euf_satask(phi)
>>> print(result)
True
"""

from sympy.assumptions.assume import global_assumptions
from sympy.assumptions.cnf import CNF, EncodedCNF
from sympy.logic.algorithms.euf_theory import EUFSolver
from sympy.logic.inference import satisfiable


def euf_satask(proposition, assumptions=True, context=global_assumptions):
    """
    Checks the satisfiability of the given proposition with assumptions using
    a SAT algorithm integrated with the EUF theory solver.

    Parameters
    ----------
    proposition : sympy.Expr
        The logical proposition (formula) to check satisfiability of.

    assumptions : sympy.Expr or bool, optional
        Extra assumptions to consider along with the proposition.

    context : sympy.assumptions.AssumptionsContext, optional
        Additional global assumptions context.

    Returns
    -------
    bool or None
        True if the proposition with assumptions is satisfiable.
        False if unsatisfiable.
        None if satisfiability unknown or undetermined under assumptions.

    Raises
    ------
    sympy.logic.algorithms.lra_theory.UnhandledInput
        If predicates not supported are encountered.

    Examples
    --------
    >>> from sympy import symbols, Eq, Function
    >>> from sympy.assumptions.euf_satask import euf_satask
    >>> f = Function('f')
    >>> a, b, x = symbols('a b x')
    >>> prop = Eq(a, b) & Eq(f(a), x)
    >>> euf_satask(prop)
    True
    """
    # Convert proposition and assumptions to CNF
    prop_cnf = CNF.from_prop(proposition)
    prop_not_cnf = CNF.from_prop(~proposition)

    assumptions_cnf = CNF.from_prop(assumptions)

    # Incorporate global assumptions context if given
    context_cnf = CNF()
    if context:
        context_cnf = context_cnf.extend(context)
    assumptions_cnf.add_from_cnf(context_cnf)

    # Encode to CNF with literal integer mapping
    coded_prop = EncodedCNF()
    coded_prop.from_cnf(prop_cnf)

    coded_not_prop = EncodedCNF()
    coded_not_prop.from_cnf(prop_not_cnf)

    coded_assumptions = EncodedCNF()
    coded_assumptions.from_cnf(assumptions_cnf)

    # Build full CNF data with assumptions and context
    cnf_true = EncodedCNF()
    cnf_true.data = coded_prop.data + coded_assumptions.data
    cnf_true.encoding = {**coded_prop.encoding, **coded_assumptions.encoding}

    cnf_false = EncodedCNF()
    cnf_false.data = coded_not_prop.data + coded_assumptions.data
    cnf_false.encoding = {**coded_not_prop.encoding, **coded_assumptions.encoding}

    # Initialize EUF theory solver on positive and negative propositions
    euf_solver_true, conflicts_true = EUFSolver.from_encoded_cnf(cnf_true)
    euf_solver_false, conflicts_false = EUFSolver.from_encoded_cnf(cnf_false)

    # Check satisfiability for prop and ~prop with assumptions
    sat_true = satisfiable(cnf_true, euf_theory=euf_solver_true) is not False
    sat_false = satisfiable(cnf_false, euf_theory=euf_solver_false) is not False

    # Interpret results consistently:
    if sat_true and not sat_false:
        return True
    elif not sat_true and sat_false:
        return False
    elif sat_true and sat_false:
        # Unknown / both models possible
        return None
    else:
        raise ValueError("Inconsistent assumptions or unexpected SAT solver results.")


from sympy import Lambda, Dummy

def preprocess_function_expression(expr):
    """
    Convert an arbitrary SymPy expression with free symbols
    into a properly curried Lambda Expression.

    Steps:
    - Extract free symbols in sorted order,
    - Create a Lambda with these symbols as variables and the original expression as body,
    - Curry the Lambda to convert multi-argument lambda into nested unary lambdas.

    Parameters
    ----------
    expr : sympy.Expr
        The expression to convert to Lambda form.

    Returns
    -------
    sympy.Lambda
        A curried Lambda representing the function.
        
    Examples
    --------
    >>> from sympy import symbols
    >>> x, y = symbols('x y')
    >>> preprocess_function_expression(x + 1)
    Lambda(x, x + 1)
    
    >>> preprocess_function_expression(x*y + y)
    Lambda(x, Lambda(y, x*y + y))
    """
    free_syms = sorted(expr.free_symbols, key=lambda s: s.name)
    if not free_syms:
        # No variables, treat as constant function with dummy argument
        dummy = Dummy('x')
        lam = Lambda(dummy, expr)
    else:
        lam = Lambda(tuple(free_syms), expr)
    # Curry the lambda to get nested 1-arg lambdas
    lam_curried = lam.curry()
    return lam_curried
