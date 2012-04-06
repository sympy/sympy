from __future__ import with_statement

import os
import tempfile
import pickle

from sympy.mpmath import *

def pickler(obj):
    fn = tempfile.mktemp()

    with open(fn, 'wb') as f:
        pickle.dump(obj, f)

    with open(fn, 'rb') as f:
        obj2 = pickle.load(f)

    os.remove(fn)

    return obj2

def test_pickle():

    obj = mpf('0.5')
    assert obj == pickler(obj)

    obj = mpc('0.5','0.2')
    assert obj == pickler(obj)
