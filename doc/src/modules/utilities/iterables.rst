=========
Iterables
=========

cartes
------

Returns the cartesian product of sequences as a generator.

Examples::
    >>> from sympy.utilities.iterables import cartes
    >>> list(cartes([1,2,3], 'ab'))
    [(1, 'a'), (1, 'b'), (2, 'a'), (2, 'b'), (3, 'a'), (3, 'b')]



variations
----------

variations(seq, n) Returns all the variations of the list of size n.

Has an optional third argument. Must be a boolean value and makes the method
return the variations with repetition if set to True, or the variations
without repetition if set to False.

Examples::
    >>> from sympy.utilities.iterables import variations
    >>> list(variations([1,2,3], 2))
    [(1, 2), (1, 3), (2, 1), (2, 3), (3, 1), (3, 2)]
    >>> list(variations([1,2,3], 2, True))
    [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (3, 3)]

Docstring
=========

.. automodule:: sympy.utilities.iterables
   :members:

