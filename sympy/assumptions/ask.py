"""Module for querying SymPy objects about assumptions."""

from sympy.assumptions.assume import (global_assumptions, Predicate,
        AppliedPredicate)
from sympy.assumptions.cnf import CNF, EncodedCNF, Literal
from sympy.core import sympify
from sympy.core.cache import cacheit
from sympy.core.kind import BooleanKind
from sympy.core.relational import Eq, Ne, Gt, Lt, Ge, Le
from sympy.logic.boolalg import (to_cnf, And, Not, Implies, Equivalent)
from sympy.logic.inference import satisfiable
from sympy.utilities.decorator import memoize_property
from sympy.utilities.exceptions import SymPyDeprecationWarning


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
        from .predicates.calculus import InfinitePredicate
        return InfinitePredicate()

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
            if args:
                facts.add(frozenset(args))
    return CNF(facts)


def ask(proposition, assumptions=True, context=global_assumptions):
    """
    Function to evaluate the proposition with assumptions.

    **Syntax**

        * ask(proposition)
            Evaluate the *proposition* in global assumption context.

        * ask(proposition, assumptions)
            Evaluate the *proposition* with respect to *assumptions* in
            global assumption context.

    This function evaluates the proposition to ``True`` or ``False`` if
    the truth value can be determined. If not, it returns ``None``.

    It should be discerned from :func:`~.refine()` which, when applied to a
    proposition, simplifies the argument to symbolic ``Boolean`` instead of
    Python built-in ``True``, ``False`` or ``None``.

    Parameters
    ==========

    proposition : any boolean expression
        Proposition which will be evaluated to boolean value. If this is
        not ``AppliedPredicate``, it will be wrapped by ``Q.is_true``.

    assumptions : any boolean expression, optional
        Local assumptions to evaluate the *proposition*.

    context : AssumptionsContext, optional
        Default assumptions to evaluate the *proposition*. By default,
        this is ``sympy.assumptions.global_assumptions`` variable.

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

    **Remarks**

        Relations in assumptions are not implemented (yet), so the following
        will not give a meaningful result.

        >>> ask(Q.positive(x), x > 0)

        It is however a work in progress.

    See Also
    ========

    sympy.assumptions.refine.refine : Simplification using assumptions.
        Proposition is not reduced to ``None`` if the truth value cannot
        be determined.
    """
    from sympy.assumptions.satask import satask

    proposition = sympify(proposition)
    assumptions = sympify(assumptions)

    if isinstance(proposition, Predicate) or proposition.kind is not BooleanKind:
        raise TypeError("proposition must be a valid logical expression")

    if isinstance(assumptions, Predicate) or assumptions.kind is not BooleanKind:
        raise TypeError("assumptions must be a valid logical expression")

    binrelpreds = {Eq: Q.eq, Ne: Q.ne, Gt: Q.gt, Lt: Q.lt, Ge: Q.ge, Le: Q.le}
    if isinstance(proposition, AppliedPredicate):
        key, args = proposition.function, proposition.arguments
    elif proposition.func in binrelpreds:
        key, args = binrelpreds[proposition.func], proposition.args
    else:
        key, args = Q.is_true, (proposition,)

    # convert local and global assumptions to CNF
    assump = CNF.from_prop(assumptions)
    assump.extend(context)

    # extract the relevant facts from assumptions with respect to args
    local_facts = _extract_all_facts(assump, args)

    known_facts_cnf = get_all_known_facts()
    known_facts_dict = get_known_facts_dict()

    # convert default facts and assumed facts to encoded CNF
    enc_cnf = EncodedCNF()
    enc_cnf.from_cnf(CNF(known_facts_cnf))
    enc_cnf.add_from_cnf(local_facts)

    # check the satisfiability of given assumptions
    if local_facts.clauses and satisfiable(enc_cnf) is False:
        raise ValueError("inconsistent assumptions %s" % assumptions)

    if local_facts.clauses:

        # quick exit if the prerequisite of proposition is not true
        # e.g. proposition = Q.odd(x), assumptions = ~Q.integer(x)
        if len(local_facts.clauses) == 1:
            cl, = local_facts.clauses
            if len(cl) == 1:
                f, = cl
                if f.is_Not and f.arg in known_facts_dict.get(key, []):
                    return False

        for clause in local_facts.clauses:
            if len(clause) == 1:
                f, = clause
                fdict = known_facts_dict.get(f.arg, None) if not f.is_Not else None
                if fdict is None:
                    pass
                elif key in fdict:
                    # quick exit if proposition is directly satisfied by assumption
                    # e.g. proposition = Q.integer(x), assumptions = Q.odd(x)
                    return True
                elif Not(key) in fdict:
                    # quick exit if proposition is directly rejected by assumption
                    # example might be proposition = Q.even(x), assumptions = Q.odd(x)
                    # but known_facts_dict does not have such information yet and
                    # such example is computed by satask.
                    return False

    # direct resolution method, no logic
    res = key(*args)._eval_ask(assumptions)
    if res is not None:
        return bool(res)
    # using satask (still costly)
    res = satask(proposition, assumptions=assumptions, context=context)
    return res


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
    class inheriting from AskHandler.

    .. deprecated:: 1.8.
        Use multipledispatch handler instead. See :obj:`~.Predicate`.

    """
    SymPyDeprecationWarning(
        feature="register_handler() function",
        useinstead="multipledispatch handler of Predicate",
        issue=20873,
        deprecated_since_version="1.8"
    ).warn()
    if isinstance(key, Predicate):
        key = key.name.name
    Qkey = getattr(Q, key, None)
    if Qkey is not None:
        Qkey.add_handler(handler)
    else:
        setattr(Q, key, Predicate(key, handlers=[handler]))


def remove_handler(key, handler):
    """
    Removes a handler from the ask system. Same syntax as register_handler

    .. deprecated:: 1.8.
        Use multipledispatch handler instead. See :obj:`~.Predicate`.

    """
    SymPyDeprecationWarning(
        feature="remove_handler() function",
        useinstead="multipledispatch handler of Predicate",
        issue=20873,
        deprecated_since_version="1.8"
    ).warn()
    if isinstance(key, Predicate):
        key = key.name.name
    getattr(Q, key).remove_handler(handler)


def single_fact_lookup(known_facts_keys, known_facts_cnf):
    # Compute the quick lookup for single facts
    mapping = {}
    for key in known_facts_keys:
        mapping[key] = {key}
        for other_key in known_facts_keys:
            if other_key != key:
                if ask_full_inference(other_key, key, known_facts_cnf):
                    mapping[key].add(other_key)
    return mapping


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
    from sympy.logic.boolalg import And
    from sympy.assumptions.cnf import Literal
    from sympy.assumptions.ask import Q

    # -{ Known facts as a set }-
    @cacheit
    def get_all_known_facts():
        return {
            %s
        }

    # -{ Known facts in Conjunctive Normal Form }-
    @cacheit
    def get_known_facts_cnf():
        return And(
            %s
        )

    # -{ Known facts in compressed sets }-
    @cacheit
    def get_known_facts_dict():
        return {
            %s
        }
    ''')
    # Compute the known facts in CNF form for logical inference
    LINE = ",\n        "
    HANG = ' '*8
    cnf = to_cnf(known_facts)
    cnf_ = CNF.to_CNF(known_facts)
    c = LINE.join([str(a) for a in cnf.args])

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
    return fact_string % (p, c, m)

@cacheit
def get_known_facts_keys():
    return [
        getattr(Q, attr)
        for attr in Q.__class__.__dict__
        if not attr.startswith('__')]

@cacheit
def get_known_facts():
    return And(
        Implies(Q.infinite, ~Q.finite),
        Implies(Q.real, Q.complex),
        Implies(Q.real, Q.hermitian),
        Equivalent(Q.extended_real, Q.real | Q.infinite),
        Equivalent(Q.even | Q.odd, Q.integer),
        Implies(Q.even, ~Q.odd),
        Implies(Q.prime, Q.integer & Q.positive & ~Q.composite),
        Implies(Q.integer, Q.rational),
        Implies(Q.rational, Q.algebraic),
        Implies(Q.algebraic, Q.complex),
        Implies(Q.algebraic, Q.finite),
        Equivalent(Q.transcendental | Q.algebraic, Q.complex & Q.finite),
        Implies(Q.transcendental, ~Q.algebraic),
        Implies(Q.transcendental, Q.finite),
        Implies(Q.imaginary, Q.complex & ~Q.real),
        Implies(Q.imaginary, Q.antihermitian),
        Implies(Q.antihermitian, ~Q.hermitian),
        Equivalent(Q.irrational | Q.rational, Q.real & Q.finite),
        Implies(Q.irrational, ~Q.rational),
        Implies(Q.zero, Q.even),

        Equivalent(Q.real, Q.negative | Q.zero | Q.positive),
        Implies(Q.zero, ~Q.negative & ~Q.positive),
        Implies(Q.negative, ~Q.positive),
        Equivalent(Q.nonnegative, Q.zero | Q.positive),
        Equivalent(Q.nonpositive, Q.zero | Q.negative),
        Equivalent(Q.nonzero, Q.negative | Q.positive),

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

from sympy.assumptions.ask_generated import (
    get_known_facts_dict, get_all_known_facts)
