"""Unification algorithm

The original version (but with some correctness issues) can be found here
http://aima.cs.berkeley.edu/python/logic.html

We extended the original algorithm with the capability to

- Handle associative-commutative unification capability.
- Designed type system to make this work with almost all Python objects by composition.

References
==========

.. [*] Artificial Intelligence: A Modern Approach by Stuart Russel and
       Peter Norvig Second edition, section 9.2, page 276
.. [*] Baader, Franz, and Tobias Nipkow. 1998. Term Rewriting and All That.
       Cambridge: Cambridge University Press. doi:10.1017/CBO9781139172752.
"""
from sympy.utilities.iterables import kbins
from abc import ABC, abstractmethod
from typing import (
    Sequence, Callable, Optional, Dict, Iterator, Union, Literal, Tuple,
    TypeVar, Any)
from itertools import permutations


T = TypeVar('T')


class Term(ABC):
    @abstractmethod
    def __str__(self) -> str:
        pass

    def __repr__(self) -> str:
        return str(self)


class Compound(Term):
    """A term with function symbol and arguments

    Attributes
    ==========

    op
        An argument to define the function head.
        Can be anything, but we recommend it to be ``str`` or something
        immutable with well-defined syntactic equality.

    args: Tuple[Term, ...]
        The arguments of the function.
        Can be a sequence of anything, but we recommend it to be a
        immutable sequence of terms.

    Examples
    ========

    How to represent functions?

    >>> from sympy.unify.core import Variable, Compound

    >>> x = Variable('x')
    >>> y = Variable('y')

    >>> Compound('f', (x,))
    f[x]
    >>> Compound('f', (x, y))
    f[x, y]

    How to represent constants?

    A 'constant' term is defined as a nullary compound term.

    >>> Compound(1, ())
    1[]
    >>> Compound(2, ())
    2[]

    For legacy reason, the following behavior is also supported

    >>> Compound('Add', (1, 2))
    Add[1, 2]

    However, the identity between the nullary compound term and its head is
    not preserved, so you need to be consistent in how constants are encoded
    into the term.

    >>> 1 == Compound(1, ())
    False
    """
    op: Any
    args: Tuple[Any, ...]

    def __init__(self, op: Any, args: Tuple[Any, ...] = ()):
        self.op = op
        self.args = tuple(args)

    def __eq__(self, other) -> bool:
        if not isinstance(other, Compound):
            return False
        if self.op != other.op:
            return False
        if len(self.args) != len(other.args):
            return False
        return all(x == y for x, y in zip(self.args, other.args))

    def __hash__(self) -> int:
        return hash((type(self), self.op, self.args))

    def __str__(self) -> str:
        return "%s[%s]" % (str(self.op), ', '.join(map(str, self.args)))


class VariableBase(Term, ABC):
    """A base class for variables"""
    pass


class CondVariable(VariableBase):
    """A variable that matches conditionally

    Explanation
    ===========

    This only has the semantics in :func:`unify`.

    Attributes
    ==========

    arg
        The name of the variable

    valid
        The function to test whether the variable should with unify with
        the term or not.
    """
    def __init__(self, arg: Any, valid: Callable[[Any], bool]):
        self.arg = arg
        self.valid = valid

    def __eq__(self, other) -> bool:
        if not isinstance(other, CondVariable):
            return False
        return self.arg == other.arg and self.valid == other.valid

    def __hash__(self) -> int:
        return hash((type(self), self.arg, self.valid))

    def __str__(self) -> str:
        return str(self.arg)


class Variable(VariableBase):
    """A variable term

    Explanation
    ===========

    It is similar as :class:`Symbol` in SymPy, however, unlike SymPy's
    :class:`Symbol` has ambiguous semantics betweeen arbitrary
    mathematical constant and functional variable, it has clear semantics
    under the substitution operation.

    Attributes
    ==========

    arg
        The name of the variable
    """
    def __init__(self, arg: Any):
        self.arg = arg

    def __eq__(self, other) -> bool:
        if not isinstance(other, Variable):
            return False
        return self.arg == other.arg

    def __hash__(self) -> int:
        return hash((type(self), self.arg))

    def __str__(self) -> str:
        return str(self.arg)


def unify(
    x: Any, y: Any,
    s: Optional[Dict[VariableBase, Any]] = None,
    **fns
) -> Iterator[Dict[VariableBase, Any]]:
    """Yields unifiers for $x$ and $y$.

    Parameters
    ==========

    x
        Any term

    y
        Any term

    s : dict, optional
        The initial substitution to be given.
        It is only needed for handling recursion internally

    is_associative : Callable[[Compound], bool], optional
        A check to whether the function should be associatively-unified.

    is_commutative : Callable[[Compound], bool], optional
        A check to whether the function should be commutatively-unified.

    Yields
    ======

    Dict[Union[Variable, CondVariable], Any]
        A unifier

    Examples
    ========

    >>> from sympy.unify.core import unify, Compound, Variable
    >>> expr = Compound("Add", ("x", "y"))
    >>> pattern = Compound("Add", ("x", Variable("a")))
    >>> next(unify(expr, pattern, {}))
    {a: 'y'}

    Notes
    =====

    This is not a standard approach of associative-commutative
    unification, but rather of list-associative-commutative unification,
    which has similar search algorithm like Unger parsing method for
    context-free grammar.

    So if you don't flatten the expressions in priori, it won't
    be able to recognize some associativity or commutativity across
    nested terms.
    """
    s = s or {}

    if x == y:
        yield s

    elif isinstance(x, VariableBase):
        yield from unify_var(x, y, s, **fns)
    elif isinstance(y, VariableBase):
        yield from unify_var(y, x, s, **fns)

    elif isinstance(x, Compound) and isinstance(y, Compound):
        is_commutative: Callable[[Compound], bool] = fns.get('is_commutative', lambda x: False)
        is_associative: Callable[[Compound], bool] = fns.get('is_associative', lambda x: False)

        if x.op != y.op:
            return

        if is_associative(x) and is_associative(y):
            a, b = (x, y) if len(x.args) < len(y.args) else (y, x)

            if is_commutative(x) and is_commutative(y):
                combs = allcombinations(a.args, b.args, 'commutative')
            else:
                combs = allcombinations(a.args, b.args, 'associative')

            for aaargs, bbargs in combs:
                aa = [unpack(Compound(a.op, arg)) for arg in aaargs]
                bb = [unpack(Compound(b.op, arg)) for arg in bbargs]
                yield from unify_multi(*zip(aa, bb), subs=s, **fns)

        elif len(x.args) == len(y.args):
            if is_commutative(x) and is_commutative(y):
                for xargs in permutations(x.args):
                    yield from unify_multi(*zip(xargs, y.args), subs=s, **fns)
            else:
                yield from unify_multi(*zip(x.args, y.args), subs=s, **fns)

    elif isinstance(x, (list, tuple)) and isinstance(y, (list, tuple)):
        if len(x) == len(y):
            yield from unify_multi(*zip(x, y), subs=s, **fns)


def unify_multi(
    *eqns: Tuple[Any, Any],
    subs: Dict[VariableBase, Any], **kwargs
) -> Iterator[Dict[VariableBase, Any]]:
    """Solve the simultaneous equation for unification"""
    if not eqns:
        yield subs
        return

    x, y = eqns[0]
    for _subs in unify(x, y, subs, **kwargs):
        yield from unify_multi(*eqns[1:], subs=_subs, **kwargs)


def unify_var(
    var: VariableBase, x: Any, s: Dict[VariableBase, Any], **fns
) -> Iterator[Dict[VariableBase, Any]]:
    if var in s:
        val = s[var]
        yield from unify(val, x, s, **fns)
    elif x in s:
        val = s[x]
        yield from unify(var, val, s, **fns)
    elif occur_check(var, x, s):
        return
    elif (isinstance(var, CondVariable) and var.valid(x)) or isinstance(var, Variable):
        s = assoc(s, var, x)
        # Cascade Substitution
        for x in s:
            s[x] = subst(s, s[x])
        yield s


def subst(s: Dict[VariableBase, Any], x: Any):
    if x in s:
        return s[x]
    if isinstance(x, Compound):
        return Compound(x.op, tuple(subst(s, arg) for arg in x.args))
    return x


def occur_check(v: VariableBase, x: Any, s: Dict[VariableBase, Any]) -> bool:
    """Checks whether the variable $v$ occurs in the term $x$"""
    if v == x:
        return True
    elif x in s:
        return occur_check(v, s[x], s)
    elif isinstance(x, Compound):
        return any(occur_check(v, arg, s) for arg in x.args)
    return False


def assoc(d, key, val):
    """ Return copy of d with key associated to val """
    d = d.copy()
    d[key] = val
    return d


def unpack(x: Any):
    if isinstance(x, Compound) and len(x.args) == 1:
        return x.args[0]
    return x


def allcombinations(
    A: Sequence[T], B: Sequence[T], ordered: str
) -> Iterator[Tuple[Tuple[Tuple[T, ...], ...], Tuple[Tuple[T, ...], ...]]]:
    """Restructure A and B to have the same number of elements.

    Parameters
    ==========

    ordered
        If ``'commutative'``, it gives all associative-commutative matches
        Otherwise, it only gives the associative matches

    Examples
    ========

    >>> from sympy.unify.core import allcombinations
    >>> for x in allcombinations((1, 2, 3), (5, 6), 'associative'): print(x)
    (((1,), (2, 3)), ((5,), (6,)))
    (((1, 2), (3,)), ((5,), (6,)))

    >>> for x in allcombinations((1, 2, 3), (5, 6), 'commutative'): print(x)
        (((1,), (2, 3)), ((5,), (6,)))
        (((1, 2), (3,)), ((5,), (6,)))
        (((1,), (3, 2)), ((5,), (6,)))
        (((1, 3), (2,)), ((5,), (6,)))
        (((2,), (1, 3)), ((5,), (6,)))
        (((2, 1), (3,)), ((5,), (6,)))
        (((2,), (3, 1)), ((5,), (6,)))
        (((2, 3), (1,)), ((5,), (6,)))
        (((3,), (1, 2)), ((5,), (6,)))
        (((3, 1), (2,)), ((5,), (6,)))
        (((3,), (2, 1)), ((5,), (6,)))
        (((3, 2), (1,)), ((5,), (6,)))
    """
    arg: Union[Literal[11], None]
    if ordered == "commutative":
        arg = 11
    else:
        arg = None

    sm, bg = (A, B) if len(A) < len(B) else (B, A)
    for part in kbins(list(range(len(bg))), len(sm), ordered=arg):
        if bg == B:
            yield tuple((a,) for a in A), partition(B, part)
        else:
            yield partition(A, part), tuple((b,) for b in B)


def partition(it, part):
    """ Partition a tuple/list into pieces defined by indices.

    Examples
    ========

    >>> from sympy.unify.core import partition
    >>> partition((10, 20, 30, 40), [[0, 1, 2], [3]])
    ((10, 20, 30), (40,))
    """
    return tuple(index(it, ind) for ind in part)


def index(it, ind):
    """Fancy indexing into an indexable iterable (tuple, list).

    Examples
    ========

    >>> from sympy.unify.core import index
    >>> index([10, 20, 30], (1, 2, 0))
    (20, 30, 10)
    """
    return tuple(it[i] for i in ind)
