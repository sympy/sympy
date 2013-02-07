from sympy import Basic
from sympy.utilities.iterables import sift

def groupby(fn, coll):
    return sift(coll, fn)

new = Basic.__new__

def is_leaf(x):
    return not isinstance(x, Basic) or x.is_Atom

def children(x):
    return x.args

def count(tup):
    return dict((k, len(v)) for k,v in groupby((lambda x: x), tup).items())
