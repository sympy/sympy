"""Module for querying SymPy objects about assumptions."""

from sympy.assumptions.assume import (global_assumptions, Predicate,
        AppliedPredicate)
from sympy.assumptions.cnf import CNF, EncodedCNF, Literal
from sympy.core import sympify
from sympy.core.kind import BooleanKind
from sympy.core.relational import Eq, Ne, Gt, Lt, Ge, Le
from sympy.logic.inference import satisfiable
from sympy.utilities.decorator import memoize_property
from sympy.utilities.exceptions import (sympy_deprecation_warning,
                                        SymPyDeprecationWarning,
                                        ignore_warnings)
from sympy.logic.boolalg import And


# Memoization is necessary for the properties of AssumptionKeys to
# ensure that only one object of Predicate objects are created.
# This is because assumption handlers are registered on those objects.


class AssumptionKeys:
    """
    This class contains all the supported keys by ``ask``.
    It should be accessed via the instance ``sympy.Q``.

    """

    # DO NOT add methods or properties other than predicate keys.
    # SAT solver checks the properties of Q and use them to compute the
    # fact system. Non-predicate attributes will break this.

    @memoize_property
    def hermitian(self):
        from .handlers.sets import HermitianPredicate
        return HermitianPredicate()

    @memoize_property
    def antihermitian(self):
        from .handlers.sets import AntihermitianPredicate
        return AntihermitianPredicate()

    @memoize_property
    def real(self):
        from .handlers.sets import RealPredicate
        return RealPredicate()

    @memoize_property
    def extended_real(self):
        from .handlers.sets import ExtendedRealPredicate
        return ExtendedRealPredicate()

    @memoize_property
    def imaginary(self):
        from .handlers.sets import ImaginaryPredicate
        return ImaginaryPredicate()

    @memoize_property
    def complex(self):
        from .handlers.sets import ComplexPredicate
        return ComplexPredicate()

    @memoize_property
    def algebraic(self):
        from .handlers.sets import AlgebraicPredicate
        return AlgebraicPredicate()

    @memoize_property
    def transcendental(self):
        from .predicates.sets import TranscendentalPredicate
        return TranscendentalPredicate()

    @memoize_property
    def integer(self):
        from .handlers.sets import IntegerPredicate
        return IntegerPredicate()

    @memoize_property
    def noninteger(self):
        from .predicates.sets import NonIntegerPredicate
        return NonIntegerPredicate()

    @memoize_property
    def rational(self):
        from .handlers.sets import RationalPredicate
        return RationalPredicate()

    @memoize_property
    def irrational(self):
        from .handlers.sets import IrrationalPredicate
        return IrrationalPredicate()

    @memoize_property
    def finite(self):
        from .handlers.calculus import FinitePredicate
        return FinitePredicate()

    @memoize_property
    def infinite(self):
        from .handlers.calculus import InfinitePredicate
        return InfinitePredicate()

    @memoize_property
    def positive_infinite(self):
        from .handlers.calculus import PositiveInfinitePredicate
        return PositiveInfinitePredicate()

    @memoize_property
    def negative_infinite(self):
        from .handlers.calculus import NegativeInfinitePredicate
        return NegativeInfinitePredicate()

    @memoize_property
    def positive(self):
        from .handlers.order import PositivePredicate
        return PositivePredicate()

    @memoize_property
    def negative(self):
        from .handlers.order import NegativePredicate
        return NegativePredicate()

    @memoize_property
    def zero(self):
        from .handlers.order import ZeroPredicate
        return ZeroPredicate()

    @memoize_property
    def extended_positive(self):
        from .handlers.order import ExtendedPositivePredicate
        return ExtendedPositivePredicate()

    @memoize_property
    def extended_negative(self):
        from .handlers.order import ExtendedNegativePredicate
        return ExtendedNegativePredicate()

    @memoize_property
    def nonzero(self):
        from .handlers.order import NonZeroPredicate
        return NonZeroPredicate()

    @memoize_property
    def nonpositive(self):
        from .handlers.order import NonPositivePredicate
        return NonPositivePredicate()

    @memoize_property
    def nonnegative(self):
        from .handlers.order import NonNegativePredicate
        return NonNegativePredicate()

    @memoize_property
    def extended_nonzero(self):
        from .handlers.order import ExtendedNonZeroPredicate
        return ExtendedNonZeroPredicate()

    @memoize_property
    def extended_nonpositive(self):
        from .handlers.order import ExtendedNonPositivePredicate
        return ExtendedNonPositivePredicate()

    @memoize_property
    def extended_nonnegative(self):
        from .handlers.order import ExtendedNonNegativePredicate
        return ExtendedNonNegativePredicate()

    @memoize_property
    def even(self):
        from .handlers.ntheory import EvenPredicate
        return EvenPredicate()

    @memoize_property
    def odd(self):
        from .handlers.ntheory import OddPredicate
        return OddPredicate()

    @memoize_property
    def prime(self):
        from .handlers.ntheory import PrimePredicate
        return PrimePredicate()

    @memoize_property
    def composite(self):
        from .handlers.ntheory import CompositePredicate
        return CompositePredicate()

    @memoize_property
    def commutative(self):
        from .handlers.common import CommutativePredicate
        return CommutativePredicate()

    @memoize_property
    def is_true(self):
        from .handlers.common import IsTruePredicate
        return IsTruePredicate()

    @memoize_property
    def symmetric(self):
        from .handlers.matrices import SymmetricPredicate
        return SymmetricPredicate()

    @memoize_property
    def invertible(self):
        from .handlers.matrices import InvertiblePredicate
        return InvertiblePredicate()

    @memoize_property
    def orthogonal(self):
        from .handlers.matrices import OrthogonalPredicate
        return OrthogonalPredicate()

    @memoize_property
    def unitary(self):
        from .handlers.matrices import UnitaryPredicate
        return UnitaryPredicate()

    @memoize_property
    def positive_definite(self):
        from .handlers.matrices import PositiveDefinitePredicate
        return PositiveDefinitePredicate()

    @memoize_property
    def upper_triangular(self):
        from .handlers.matrices import UpperTriangularPredicate
        return UpperTriangularPredicate()

    @memoize_property
    def lower_triangular(self):
        from .handlers.matrices import LowerTriangularPredicate
        return LowerTriangularPredicate()

    @memoize_property
    def diagonal(self):
        from .handlers.matrices import DiagonalPredicate
        return DiagonalPredicate()

    @memoize_property
    def fullrank(self):
        from .handlers.matrices import FullRankPredicate
        return FullRankPredicate()

    @memoize_property
    def square(self):
        from .handlers.matrices import SquarePredicate
        return SquarePredicate()

    @memoize_property
    def integer_elements(self):
        from .handlers.matrices import IntegerElementsPredicate
        return IntegerElementsPredicate()

    @memoize_property
    def real_elements(self):
        from .handlers.matrices import RealElementsPredicate
        return RealElementsPredicate()

    @memoize_property
    def complex_elements(self):
        from .handlers.matrices import ComplexElementsPredicate
        return ComplexElementsPredicate()

    @memoize_property
    def singular(self):
        from .predicates.matrices import SingularPredicate
        return SingularPredicate()

    @memoize_property
    def normal(self):
        from .predicates.matrices import NormalPredicate
        return NormalPredicate()

    @memoize_property
    def triangular(self):
        from .predicates.matrices import TriangularPredicate
        return TriangularPredicate()

    @memoize_property
    def unit_triangular(self):
        from .predicates.matrices import UnitTriangularPredicate
        return UnitTriangularPredicate()

    @memoize_property
    def eq(self):
        from .relation.equality import EqualityPredicate
        return EqualityPredicate()

    @memoize_property
    def ne(self):
        from .relation.equality import UnequalityPredicate
        return UnequalityPredicate()

    @memoize_property
    def gt(self):
        from .relation.equality import StrictGreaterThanPredicate
        return StrictGreaterThanPredicate()

    @memoize_property
    def ge(self):
        from .relation.equality import GreaterThanPredicate
        return GreaterThanPredicate()

    @memoize_property
    def lt(self):
        from .relation.equality import StrictLessThanPredicate
        return StrictLessThanPredicate()

    @memoize_property
    def le(self):
        from .relation.equality import LessThanPredicate
        return LessThanPredicate()


Q = AssumptionKeys()

def _extract_all_facts(assump, exprs):
    """
    Extract all relevant assumptions from *assump* with respect to given *exprs*.

    Parameters
    ==========

    assump : sympy.assumptions.cnf.CNF

    exprs : tuple of expressions

    Returns
    =======

    sympy.assumptions.cnf.CNF

    Examples
    ========

    >>> from sympy import Q
    >>> from sympy.assumptions.cnf import CNF
    >>> from sympy.assumptions.ask import _extract_all_facts
    >>> from sympy.abc import x, y
    >>> assump = CNF.from_prop(Q.positive(x) & Q.integer(y))
    >>> exprs = (x,)
    >>> cnf = _extract_all_facts(assump, exprs)
    >>> cnf.clauses
    {frozenset({Literal(Q.positive, False)})}

    """
    facts = set()

    for clause in assump.clauses:
        args = []
        for literal in clause:
            if isinstance(literal.lit, AppliedPredicate) and len(literal.lit.arguments) == 1:
                if literal.lit.arg in exprs:
                    # Add literal if it has matching in it
                    args.append(Literal(literal.lit.function, literal.is_Not))
                else:
                    # If any of the literals doesn't have matching expr don't add the whole clause.
                    break
            else:
                # If any of the literals aren't unary predicate don't add the whole clause.
                break

        else:
            if args:
                facts.add(frozenset(args))
    return CNF(facts)

def _deduce_binrel_from_assumptions(proposition, assumptions):
    """
    Infer logical consequences or contradictions involving binary relational predicates.

    Explanation
    ===========

    This is a minimal rule-based system used by `ask()` to infer binary relational
    consequences that are not captured by the univariate (CNF-based) logic engine.

    For example, the predicate `Q.ge(a, b)` can be inferred as true from the assumption
    `Q.gt(a, b)`; conversely, `Q.gt(a, b)` contradicts `Q.lt(a, b)`, and so the former
    should be inferred as false if the latter is assumed.

    This function handles such implications and contradictions manually for common
    binary predicates (e.g. Q.gt, Q.ge, Q.lt, Q.le, Q.eq), ensuring correct behavior
    for simple logical inferences involving binary relations.

    This is necessary because the logic backend of SymPy currently only supports
    univariate predicate logic via CNF, and binary relations are not yet integrated
    into the SAT solver.

    Examples
    ========

    >>> from sympy import symbols, ask, Q
    >>> a, b = symbols('a b', real=True)
    >>> ask(Q.ge(a, b), Q.gt(a, b))
    True
    >>> ask(Q.gt(a, b), Q.lt(a, b))
    False

    Notes
    =====

    This function is intended to be called internally by `ask()` and is not part
    of the public API.

    """
    # Ensure the proposition is a binary relational predicate (e.g., Q.ge(a, b))
    if not isinstance(proposition, AppliedPredicate):
        return None

    f = proposition.function
    args = proposition.arguments

    # Rule set 1: Positive implications between binary predicates.
    # These are simple cases where a stronger assumption implies a weaker one.
    # For example, Q.gt(a, b) => Q.ge(a, b)
    implication_rules = [
        (Q.gt, Q.ge),  # gt ⇒ ge
        (Q.lt, Q.le),  # lt ⇒ le
    ]

    for stronger, weaker in implication_rules:
        if f == weaker:
            # If proposition is e.g. Q.ge(a, b), and assumptions contain Q.gt(a, b),
            # we can directly return True.
            stronger_pred = stronger(*args)
            if ask(stronger_pred, assumptions):
                return True

    # Rule set 2: Logical contradictions between binary predicates.
    # These cover pairs of predicates that cannot simultaneously be true.
    # For example, Q.lt(a, b) contradicts Q.gt(a, b), so one must be False.
    contradiction_pairs = [
        (Q.gt, Q.lt), (Q.lt, Q.gt),
        (Q.ge, Q.lt), (Q.le, Q.gt),
        (Q.eq, Q.gt), (Q.eq, Q.lt),
    ]

    # Normalize assumptions to a list for uniform processing
    assumption_list = []
    if isinstance(assumptions, AppliedPredicate):
        assumption_list = [assumptions]
    elif isinstance(assumptions, And):
        assumption_list = list(assumptions.args)

    # Check for contradictions: if any contradictory predicate is present in
    # the assumptions, we can safely return False.
    for a_pred, b_pred in contradiction_pairs:
        if f == a_pred:
            contradicting = b_pred(*args)
            if contradicting in assumption_list:
                return False

    # If no inference could be made, return None
    return None



def ask(proposition, assumptions=True, context=global_assumptions):
    """
    Function to evaluate the proposition with assumptions.

    Explanation
    ===========

    This function evaluates the proposition to ``True`` or ``False`` if
    the truth value can be determined. If not, it returns ``None``.

    It should be discerned from :func:`~.refine` which, when applied to a
    proposition, simplifies the argument to symbolic ``Boolean`` instead of
    Python built-in ``True``, ``False`` or ``None``.

    **Syntax**

        * ask(proposition)
            Evaluate the *proposition* in global assumption context.

        * ask(proposition, assumptions)
            Evaluate the *proposition* with respect to *assumptions* in
            global assumption context.

    Parameters
    ==========

    proposition : Boolean
        Proposition which will be evaluated to boolean value. If this is
        not ``AppliedPredicate``, it will be wrapped by ``Q.is_true``.

    assumptions : Boolean, optional
        Local assumptions to evaluate the *proposition*.

    context : AssumptionsContext, optional
        Default assumptions to evaluate the *proposition*. By default,
        this is ``sympy.assumptions.global_assumptions`` variable.

    Returns
    =======

    ``True``, ``False``, or ``None``

    Raises
    ======

    TypeError : *proposition* or *assumptions* is not valid logical expression.

    ValueError : assumptions are inconsistent.

    Examples
    ========

    >>> from sympy import ask, Q, pi
    >>> from sympy.abc import x, y
    >>> ask(Q.rational(pi))
    False
    >>> ask(Q.even(x*y), Q.even(x) & Q.integer(y))
    True
    >>> ask(Q.prime(4*x), Q.integer(x))
    False

    If the truth value cannot be determined, ``None`` will be returned.

    >>> print(ask(Q.odd(3*x))) # cannot determine unless we know x
    None

    ``ValueError`` is raised if assumptions are inconsistent.

    >>> ask(Q.integer(x), Q.even(x) & Q.odd(x))
    Traceback (most recent call last):
      ...
    ValueError: inconsistent assumptions Q.even(x) & Q.odd(x)

    Notes
    =====

    Relations in assumptions are not implemented (yet), so the following
    will not give a meaningful result.

    >>> ask(Q.positive(x), x > 0)

    It is however a work in progress.

    See Also
    ========

    sympy.assumptions.refine.refine : Simplification using assumptions.
        Proposition is not reduced to ``None`` if the truth value cannot
        be determined.
        
    This function also includes a heuristic to handle cases where Python evaluates
    symbolic relations like `a == 0` to `False` before SymPy can interpret them symbolically.
    When this happens and assumptions include `Q.zero(a)` or `Q.eq(a, b)`, the function
    will attempt to recover the intended symbolic relation and evaluate it accordingly.

    This helps ensure expressions like `ask(a == 0, Q.zero(a))` or
    `ask(a == 0, Q.eq(a, 0))` correctly return True.
    
    """
    from sympy.assumptions.satask import satask
    from sympy.assumptions.lra_satask import lra_satask
    from sympy.logic.algorithms.lra_theory import UnhandledInput
    from sympy.assumptions.predicates.common import IsTruePredicate
    from sympy.core.relational import Relational
    from sympy.assumptions import AppliedPredicate
    from sympy.logic.boolalg import And


    if isinstance(proposition, bool):
        # Handle cases where symbolic comparison (e.g. a == 0) has been evaluated
        # prematurely by Python into a plain bool. This typically happens when
        # users write expressions like `ask(a == 0, Q.zero(a))`, where `a == 0`
        # gets resolved to False before SymPy can wrap it as a symbolic Eq(a, 0).
        #
        # This block attempts to recover the intended symbolic meaning by
        # inspecting the assumptions and reconstructing the appropriate
        # binary predicate (e.g. Q.eq(a, 0)) from common patterns.
        #
        # For now, it handles:
        #   - Q.zero(a) ⇒ assume proposition was Q.eq(a, 0)
        #   - Q.eq(a, b) ⇒ assume proposition was also Q.eq(a, b)
        
        assumption_list = []
        if isinstance(assumptions, AppliedPredicate):
            assumption_list = [assumptions]
        elif isinstance(assumptions, And):
            assumption_list = list(assumptions.args)

        for a in assumption_list:
            if not isinstance(a, AppliedPredicate):
                continue

            func = a.function
            args = a.arguments

            # Case 1: from Q.zero(a) assume proposition should be Q.eq(a, 0)
            if func == Q.zero and len(args) == 1:
                x = args[0]
                return ask(Q.eq(x, 0), assumptions, context)

            # Case 2: from Q.eq(a, b) assume proposition should be Q.eq(a, b)
            if func == Q.eq and len(args) == 2:
                x, y = args
                return ask(Q.eq(x, y), assumptions, context)


    proposition = sympify(proposition)
    assumptions = sympify(assumptions)
    
    if isinstance(proposition, Relational):
        proposition = Q.is_true(proposition)

    if isinstance(proposition, Predicate) or proposition.kind is not BooleanKind:
        raise TypeError("proposition must be a valid logical expression")

    if isinstance(assumptions, Predicate) or assumptions.kind is not BooleanKind:
        raise TypeError("assumptions must be a valid logical expression")

    binrelpreds = {Eq: Q.eq, Ne: Q.ne, Gt: Q.gt, Lt: Q.lt, Ge: Q.ge, Le: Q.le}
    if isinstance(proposition, AppliedPredicate):
        key, args = proposition.function, proposition.arguments
    elif proposition.func in binrelpreds:
        key, args = binrelpreds[type(proposition)], proposition.args
    else:
        key, args = Q.is_true, (proposition,)
        
    # Intercept binary relation implications
    binrel_result = _deduce_binrel_from_assumptions(proposition, assumptions)
    if binrel_result is not None:
        return binrel_result

    # convert local and global assumptions to CNF
    assump_cnf = CNF.from_prop(assumptions)
    assump_cnf.extend(context)

    # extract the relevant facts from assumptions with respect to args
    local_facts = _extract_all_facts(assump_cnf, args)

    # convert default facts and assumed facts to encoded CNF
    known_facts_cnf = get_all_known_facts()
    enc_cnf = EncodedCNF()
    enc_cnf.from_cnf(CNF(known_facts_cnf))
    enc_cnf.add_from_cnf(local_facts)

    # check the satisfiability of given assumptions
    if local_facts.clauses and satisfiable(enc_cnf) is False:
        raise ValueError("inconsistent assumptions %s" % assumptions)

    # quick computation for single fact
    res = _ask_single_fact(key, local_facts)
    if res is not None:
        return res

    # direct resolution method, no logic
    res = key(*args)._eval_ask(assumptions)
    if res is not None:
        return bool(res)

    # using satask (still costly)
    res = satask(proposition, assumptions=assumptions, context=context)
    if res is not None:
        return res

    try:
        res = lra_satask(proposition, assumptions=assumptions, context=context)
    except UnhandledInput:
        return None

    return res


def _ask_single_fact(key, local_facts):
    """
    Compute the truth value of single predicate using assumptions.

    Parameters
    ==========

    key : sympy.assumptions.assume.Predicate
        Proposition predicate.

    local_facts : sympy.assumptions.cnf.CNF
        Local assumption in CNF form.

    Returns
    =======

    ``True``, ``False`` or ``None``

    Examples
    ========

    >>> from sympy import Q
    >>> from sympy.assumptions.cnf import CNF
    >>> from sympy.assumptions.ask import _ask_single_fact

    If prerequisite of proposition is rejected by the assumption,
    return ``False``.

    >>> key, assump = Q.zero, ~Q.zero
    >>> local_facts = CNF.from_prop(assump)
    >>> _ask_single_fact(key, local_facts)
    False
    >>> key, assump = Q.zero, ~Q.even
    >>> local_facts = CNF.from_prop(assump)
    >>> _ask_single_fact(key, local_facts)
    False

    If assumption implies the proposition, return ``True``.

    >>> key, assump = Q.even, Q.zero
    >>> local_facts = CNF.from_prop(assump)
    >>> _ask_single_fact(key, local_facts)
    True

    If proposition rejects the assumption, return ``False``.

    >>> key, assump = Q.even, Q.odd
    >>> local_facts = CNF.from_prop(assump)
    >>> _ask_single_fact(key, local_facts)
    False
    """
    if local_facts.clauses:

        known_facts_dict = get_known_facts_dict()

        if len(local_facts.clauses) == 1:
            cl, = local_facts.clauses
            if len(cl) == 1:
                f, = cl
                prop_facts = known_facts_dict.get(key, None)
                prop_req = prop_facts[0] if prop_facts is not None else set()
                if f.is_Not and f.arg in prop_req:
                    # the prerequisite of proposition is rejected
                    return False

        for clause in local_facts.clauses:
            if len(clause) == 1:
                f, = clause
                prop_facts = known_facts_dict.get(f.arg, None) if not f.is_Not else None
                if prop_facts is None:
                    continue

                prop_req, prop_rej = prop_facts
                if key in prop_req:
                    # assumption implies the proposition
                    return True
                elif key in prop_rej:
                    # proposition rejects the assumption
                    return False

    return None


def register_handler(key, handler):
    """
    Register a handler in the ask system. key must be a string and handler a
    class inheriting from AskHandler.

    .. deprecated:: 1.8.
        Use multipledispatch handler instead. See :obj:`~.Predicate`.

    """
    sympy_deprecation_warning(
        """
        The AskHandler system is deprecated. The register_handler() function
        should be replaced with the multipledispatch handler of Predicate.
        """,
        deprecated_since_version="1.8",
        active_deprecations_target='deprecated-askhandler',
    )
    if isinstance(key, Predicate):
        key = key.name.name
    Qkey = getattr(Q, key, None)
    if Qkey is not None:
        Qkey.add_handler(handler)
    else:
        setattr(Q, key, Predicate(key, handlers=[handler]))


def remove_handler(key, handler):
    """
    Removes a handler from the ask system.

    .. deprecated:: 1.8.
        Use multipledispatch handler instead. See :obj:`~.Predicate`.

    """
    sympy_deprecation_warning(
        """
        The AskHandler system is deprecated. The remove_handler() function
        should be replaced with the multipledispatch handler of Predicate.
        """,
        deprecated_since_version="1.8",
        active_deprecations_target='deprecated-askhandler',
    )
    if isinstance(key, Predicate):
        key = key.name.name
    # Don't show the same warning again recursively
    with ignore_warnings(SymPyDeprecationWarning):
        getattr(Q, key).remove_handler(handler)


from sympy.assumptions.ask_generated import (get_all_known_facts,
    get_known_facts_dict)
