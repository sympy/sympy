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

class RuleDB(object):
    def __init__(self):
        self.coll = list()

    def insert(self, tup, val):
        self.coll.append((count(tup), val))

    def query(self, indata):
        inkey = count(indata)
        for collkey, collval in self.coll:
            if all(inkey.get(k, 0) >= v for k, v in collkey.items()):
                yield collval
