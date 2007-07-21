
def all(iterable):
    """Return True if all elements are set to True. This
       function does not support predicates explicitely,
       but this behaviour can be simulated easily using
       list comprehension.

       >>> from sympy.modules.utilities import all

       >>> all( [True, True, True] )
       True
       >>> all( [True, False, True] )
       False
       >>> all( [ x % 2 == 0 for x in [2, 6, 8] ] )
       True
       >>> all( [ x % 2 == 0 for x in [2, 6, 7] ] )
       False

       NOTE: Starting from Python 2.5 this a built-in.
    """
    for item in iterable:
        if not item:
            return False
    return True

def any(iterable):
    """Return True if at least one element is set to True.
       This function does not support predicates explicitely,
       but this behaviour can be simulated easily using
       list comprehension.

       >>> from sympy.modules.utilities import any
       >>> any( [False, False, False] )
       False
       >>> any( [False, True, False] )
       True
       >>> any( [ x % 2 == 1 for x in [2, 6, 8] ] )
       False
       >>> any( [ x % 2 == 1 for x in [2, 6, 7] ] )
       True

       NOTE: Starting from Python 2.5 this a built-in.
    """
    for item in iterable:
        if item:
            return True
    return False

def make_list(expr, kind):
    """Returns a list of elements taken from specified expresion
       when it is of sequence type (Add or Mul) or singleton list
       otherwise (Rational, Pow etc.).

       >>> from sympy import *
       >>> x, y = map(Symbol, 'xy')
       >>> make_list(x*y, Mul)
       [x, y]
       >>> make_list(x*y, Add)
       [x * y]
       >>> make_list(x*y + y, Add)
       [y, x * y]

    """
    if isinstance(expr, kind):
        return list(expr[:])
    else:
        return [expr]

def flatten(iterable):
    """Recursively denest iterable containers.

       >>> flatten([1, 2, 3])
       [1, 2, 3]
       >>> flatten([1, 2, [3]])
       [1, 2, 3]
       >>> flatten([1, [2, 3], [4, 5]])
       [1, 2, 3, 4, 5]
    """
    result = []

    for item in iterable:
        if isinstance(item, list):
            result.extend(flatten(item))
        else:
            result.append(item)

    return result
