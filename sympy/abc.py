"""
This module exports all latin and greek letters as Symbols, so you can
conveniently do

    >>> from sympy.abc import x, y

instead of the slightly more clunky-looking

    >>> from sympy import symbols
    >>> x, y = symbols('x y')

Caveats
=======

1. As of the time of writing this, the names ``C``, ``O``, ``S``, ``I``, ``N``,
``E``, and ``Q`` are colliding with names defined in SymPy. If you import them
both from ``sympy.abc`` and ``sympy``, the second import will "win".
I.e. you can do ``from sympy.abc import *; from sympy import *`` to get all
``sympy.abc`` symbols and still have ``C`` with the usual definition from SymPy.
However, this approach is recommended only for interactive sessions and
short-lived code, because future versions of ``from sympy import *`` may import
more single-letter names than today.

2. This module does not define symbol names on demand, i.e.
```from sympy.abc import foo``` will be reported as an error because
``sympy.abc`` does not generate the name ``foo``. To get a symbol named `'foo'`,
you still need to use ``Symbol('foo')`` or ``symbols('foo')``.
You can freely mix usage of ``sympy.abc`` and ``Symbol``/``symbols``, though if
the code is supposed to live longer than the current session, it's usually
easier on human readers to stick with one and only one way to get the required
symbols.

3. Using ``from sympy.abc import x, y`` will confuse static analysis tools
because they won't find the definition of ``x`` or ``y`` in the source of
``sympy/abc.py``.
This is more of a problem for larger amounts of ``abc``-using code, such as
SymPy itself.
"""

from __future__ import print_function, division

import string

from .core import Symbol
from .core.alphabets import greeks
from .core.compatibility import exec_

_latin = list(string.ascii_letters)
# COSINEQ should not be imported as they clash; gamma, pi and zeta clash, too
_greek = list(greeks) # make a copy, so we can mutate it
# Note: We import lamda since lambda is a reserved keyword in Python
_greek.remove("lambda")
_greek.append("lamda")

for _s in _latin + _greek:
    exec_("%s = Symbol('%s')" % (_s, _s))

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
    exec_('from sympy import *', ns)
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
