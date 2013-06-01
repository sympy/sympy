from core import Symbol

_latin = list('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')
# COSINEQ should not be imported as they clash; gamma, pi and zeta clash, too
_greek = 'alpha beta gamma delta epsilon zeta eta theta iota kappa lamda '\
    'mu nu xi omicron pi rho sigma tau upsilon phi chi psi omega'.split(' ')
# Note: We import lamda since lambda is a reserved keyword in Python

for _s in _latin + _greek:
    exec "%s = Symbol('%s')" % (_s, _s)

def clashing():
    """Return the clashing-symbols dictionaries.

    ``clash1`` defines all the single letter variables that clash with
    SymPy objects; ``clash2`` defines the multi-letter clashing symbols;
    and ``clash`` is the union of both. These can be passed for ``locals``
    during sympification if one desires Symbols rather than the non-Symbol
    objects for those names.

    Examples
    ========

    >>> from sympy import S
    >>> from sympy.abc import _clash1, _clash2, _clash
    >>> S("Q & C", locals=_clash1)
    And(C, Q)
    >>> S('pi(x)', locals=_clash2)
    pi(x)
    >>> S('pi(C, Q)', locals=_clash)
    pi(C, Q)

    Note: if changes are made to the docstring examples they can only
    be tested after removing "clashing" from the list of deleted items
    at the bottom of this file which removes this function from the
    namespace.
    """

    ns = {}
    exec 'from sympy import *' in ns
    clash1 = {}
    clash2 = {}
    while ns:
        k, _ = ns.popitem()
        if k in _greek:
            clash2[k] = Symbol(k)
            _greek.remove(k)
        elif k in _latin:
            clash1[k] = Symbol(k)
            _latin.remove(k)
    clash = {}
    clash.update(clash1)
    clash.update(clash2)
    return clash1, clash2, clash

_clash1, _clash2, _clash = clashing()

del _latin, _greek, _s, clashing, Symbol
