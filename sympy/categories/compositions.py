"""
Makes composition of any object possible in a mathematical way.
f( g( x ) ) = (f*g)(x)
f(g) = f*g

Illustration:
============
Let, f,g,h,...n be callable objects
then (f*g*h*...*n) will be a `Compose` object
and the corresponding native object would `(f*g*h*...*n).func`.
This notation makes it easy to do category theory in its natural state in python.


sympy example:
==============
>>> from sympy import Symbol, sin, cos, diff
>>> x= Symbol('x')
>>> f = Compose(lambda x:sin(x))
>>> g = Compose(lambda x:cos(x))
>>> y = (f*g*f)(x)
>>> print(diff(y, x))
-sin(sin(x))*cos(x)*cos(cos(sin(x)))
"""


class Compose(object):
    """
    Makes function composable
    >>> @Compose
    ... def f(x):
    ...     return x
    ...
    >>> g = Compose(lambda x: x)
    >>> print((f * g)(2))
    2
    >>> from sympy import Symbol, sin, cos, diff
    >>> x= Symbol('x')
    >>> f = Compose(lambda x:sin(x))
    >>> g = Compose(lambda x:cos(x))
    >>> y = (f*g*f)(x)
    >>> print(diff(y, x))
    -sin(sin(x))*cos(x)*cos(cos(sin(x)))
    """

    def __init__(self, func):
        self.func = func

    def __call__(self, x):
        return self.func(x)

    def __mul__(self, neighbour):
        return Compose(lambda x: self.func(neighbour.func(x)))


def tests():
    # Syntax 1
    @Compose
    def f(x):
        return x

    # Syntax 2
    g = Compose(lambda x: x)

    print((f * g)(2))

    @Compose
    def log(text):
        return "Log " + text

    @Compose
    def warn(text):
        return " Warn " + text

    print((warn * log)("this"))


if __name__ == "__main__":
    tests()
    import doctest

    doctest.testmod()
