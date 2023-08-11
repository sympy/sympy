from sympy.core.function import AppliedUndef
from sympy.sets.fancysets import ImageSet
from sympy.sets.sets import FiniteSet
from sympy.tensor.indexed import Indexed


def _get_free_symbols(exprs):
    """Returns the free symbols of a symbolic expression.

    If the expression contains any of these elements, assume that they are
    the "free symbols" of the expression:

    * indexed objects
    * applied undefined function (useful for sympy.physics.mechanics module)
    """
    if not isinstance(exprs, (list, tuple, set)):
        exprs = [exprs]
    if all(callable(e) for e in exprs):
        return set()

    free = set().union(*[e.atoms(Indexed) for e in exprs])
    free = free.union(*[e.atoms(AppliedUndef) for e in exprs])
    if len(free) > 0:
        return free
    return set().union(*[e.free_symbols for e in exprs])


def extract_solution(set_sol, n=10):
    """Extract numerical solutions from a set solution (computed by solveset,
    linsolve, nonlinsolve). Often, it is not trivial do get something useful
    out of them.

    Parameters
    ==========

    n : int, optional
        In order to replace ImageSet with FiniteSet, an iterator is created
        for each ImageSet contained in `set_sol`, starting from 0 up to `n`.
        Default value: 10.
    """
    images = set_sol.find(ImageSet)
    for im in images:
        it = iter(im)
        s = FiniteSet(*[next(it) for n in range(0, n)])
        set_sol = set_sol.subs(im, s)
    return set_sol
