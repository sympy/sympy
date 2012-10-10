# Generic rules for SymPy.
# This file assumes knowledge of Basic and little else
from sympy import Basic

def rm_id(isid):
    """ Create a rule to remove identities

    isid - fn :: x -> Bool  --- whether or not this element is an identity

    >>> from sympy.rr import rm_id, typed
    >>> from sympy import Basic
    >>> remove_zeros = rm_id(lambda x: x==0)
    >>> remove_zeros(Basic(1, 0, 2))
    Basic(1, 2)
    >>> remove_zeros(Basic(0, 0))
    Basic(0)
    """
    def ident_remove(expr):
        """ Remove identities """
        ids = map(isid, expr.args)
        if sum(ids) == 0:           # No identities. Common case
            return expr
        elif sum(ids) != len(ids):  # there is at least one non-identity
            return Basic.__new__(expr.__class__,
                              *[arg for arg, x in zip(expr.args, ids) if not x])
        else:
            first_id = (arg for arg, x in zip(expr.args, ids) if x).next()
            return Basic.__new__(expr.__class__, first_id)

    return ident_remove

def frequencies(coll):
    """ Count the number of times an element occurs in a collection

    >>> from sympy.rr.rl import frequencies
    >>> frequencies("aabbc")
    {'a': 2, 'b': 2, 'c': 1}
    """
    counts = {}
    for elem in coll:
        counts[elem] = counts.get(elem, 0) + 1
    return counts

def glom(mkglom):
    """ Create a rule to conglomerate identical args

    >>> from sympy.rr import glom
    >>> from sympy import Basic
    >>> rl = glom(lambda num, arg: num * arg)
    >>> rl(Basic(1, 1, 3))
    Basic(2, 3)
    """
    def conglomerate(expr):
        """ Conglomerate together identical args x + x -> 2x """
        freqs = frequencies(expr.args)
        return Basic.__new__(type(expr), *[arg if freqs[arg] == 1
                                               else mkglom(freqs[arg], arg)
                                               for arg in freqs])
    return conglomerate

def unpack(expr):
    """ Rule to unpack singleton args

    >>> from sympy.rr import unpack
    >>> from sympy import Basic
    >>> unpack(Basic(2))
    2
    """
    if len(expr.args) == 1:
        return expr.args[0]
    else:
        return expr

def flatten(expr):
    """ Flatten T(a, b, T(c, d), T2(e)) to T(a, b, c, d, T2(e)) """
    cls = expr.__class__
    args = []
    for arg in expr.args:
        if arg.__class__ == cls:
            args.extend(arg.args)
        else:
            args.append(arg)
    return Basic.__new__(expr.__class__, *args)

def sort(key):
    """ Create a rule to sort by a key function

    >>> from sympy.rr import sort
    >>> from sympy import Basic
    >>> sort_rl = sort(str)
    >>> sort_rl(Basic(3, 1, 2))
    Basic(1, 2, 3)
    """

    def sort_rl(expr):
        return Basic.__new__(expr.__class__, *sorted(expr.args, key=key))
    return sort_rl

def distribute(A, B):
    """ Turns an A containing Bs into a B of As

    where A, B are container types

    >>> from sympy.rr import distribute
    >>> from sympy import Add, Mul, symbols
    >>> x, y = symbols('x,y')
    >>> dist = distribute(Mul, Add)
    >>> expr = Mul(2, x+y, evaluate=False)
    >>> expr
    2*(x + y)
    >>> dist(expr)
    2*x + 2*y
    """

    def distribute_rl(expr):
        if not any(isinstance(arg, B) for arg in expr.args):
            return expr

        for i, arg in enumerate(expr.args):
            if isinstance(arg, B): break

        first, b, tail = expr.args[:i], expr.args[i], expr.args[i+1:]
        return B(*[A(*(first + (arg,) + tail)) for arg in b.args])
    return distribute_rl
