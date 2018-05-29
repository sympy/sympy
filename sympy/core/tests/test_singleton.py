from sympy.core.basic import Basic
from sympy.core.numbers import Rational
from sympy.core.singleton import S, Singleton, SingletonRegistry

from sympy.core.compatibility import with_metaclass, exec_

def test_Singleton():
    global instantiated
    instantiated = 0

    class MySingleton(with_metaclass(Singleton, Basic)):
        def __new__(cls):
            global instantiated
            instantiated += 1
            return Basic.__new__(cls)

    assert instantiated == 0
    MySingleton() # force instantiation
    assert instantiated == 1
    assert MySingleton() is not Basic()
    assert MySingleton() is MySingleton()
    assert S.MySingleton is MySingleton()
    assert instantiated == 1

    class MySingleton_sub(MySingleton):
        pass
    assert instantiated == 1
    MySingleton_sub()
    assert instantiated == 2
    assert MySingleton_sub() is not MySingleton()
    assert MySingleton_sub() is MySingleton_sub()

def test_names_in_namespace():
    # Every singleton name should be accessible from the 'from sympy import *'
    # namespace in addition to the S object. However, it does not need to be
    # by the same name (e.g., oo instead of S.Infinity).
    d = {}
    exec_('from sympy import *', d)

    for name in dir(S):
        if name.startswith('_'):
            continue
        if hasattr(SingletonRegistry, name):
            continue
        if isinstance(getattr(S, name), Rational):
            continue
        if name == 'MySingleton':
            # From the test above
            continue
        if name == 'NegativeInfinity':
            # Accessible by -oo
            continue

        # Use is here because of complications with ==
        assert any(getattr(S, name) is i for i in d.values()), name
