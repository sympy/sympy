"""sympify -- convert objects SymPy internal format"""

from inspect import getmro

from core import all_classes as sympy_classes
from sympy.core.compatibility import iterable


class SympifyError(ValueError):
    def __init__(self, expr, base_exc=None):
        self.expr = expr
        self.base_exc = base_exc
    def __str__(self):
        if self.base_exc is None:
            return "SympifyError: %r" % (self.expr,)

        return ("Sympify of expression '%s' failed, because of exception being "
            "raised:\n%s: %s" % (self.expr, self.base_exc.__class__.__name__,
            str(self.base_exc)))

converter = {} #See sympify docstring.

def sympify(a, locals=None, convert_xor=True, strict=False, rational=False):
    """
    Converts an arbitrary expression to a type that can be used inside sympy.

    For example, it will convert python ints into instance of sympy.Rational,
    floats into instances of sympy.Float, etc. It is also able to coerce symbolic
    expressions which inherit from Basic. This can be useful in cooperation
    with SAGE.

    It currently accepts as arguments:
       - any object defined in sympy (except matrices [TODO])
       - standard numeric python types: int, long, float, Decimal
       - strings (like "0.09" or "2e-19")
       - booleans, including ``None`` (will leave them unchanged)
       - lists, sets or tuples containing any of the above

    If the argument is already a type that sympy understands, it will do
    nothing but return that value. This can be used at the beginning of a
    function to ensure you are working with the correct type.

    >>> from sympy import sympify

    >>> sympify(2).is_integer
    True
    >>> sympify(2).is_real
    True

    >>> sympify(2.0).is_real
    True
    >>> sympify("2.0").is_real
    True
    >>> sympify("2e-45").is_real
    True

    If the expression could not be converted, a SympifyError is raised.

    >>> sympify("x***2")
    Traceback (most recent call last):
    ...
    SympifyError: SympifyError: "could not parse u'x***2'"


    If the option ``strict`` is set to ``True``, only the types for which an
    explicit conversion has been defined are converted. In the other
    cases, a SympifyError is raised.

    >>> sympify(True)
    True
    >>> sympify(True, strict=True)
    Traceback (most recent call last):
    ...
    SympifyError: SympifyError: True

    To extend ``sympify`` to convert custom objects (not derived from ``Basic``),
    just define a ``_sympy_`` method to your class. You can do that even to
    classes that you do not own by subclassing or adding the method at runtime.

    >>> from sympy import Matrix
    >>> class MyList1(object):
    ...     def __iter__(self):
    ...         yield 1
    ...         yield 2
    ...         raise StopIteration
    ...     def __getitem__(self, i): return list(self)[i]
    ...     def _sympy_(self): return Matrix(self)
    >>> sympify(MyList1())
    [1]
    [2]

    If you do not have control over the class definition you could also use the
    ``converter`` global dictionary. The key is the class and the value is a
    function that takes a single argument and returns the desired SymPy
    object, e.g. ``converter[MyList] = lambda x: Matrix(x)``.

    >>> class MyList2(object):   # XXX Do not do this if you control the class!
    ...     def __iter__(self):  #     Use _sympy_!
    ...         yield 1
    ...         yield 2
    ...         raise StopIteration
    ...     def __getitem__(self, i): return list(self)[i]
    >>> from sympy.core.sympify import converter
    >>> converter[MyList2] = lambda x: Matrix(x)
    >>> sympify(MyList2())
    [1]
    [2]

    """
    try:
        cls = a.__class__
    except AttributeError:  #a is probably an old-style class object
        cls = type(a)
    if cls in sympy_classes:
        return a
    if cls in (bool, type(None)):
        if strict:
            raise SympifyError(a)
        else:
            return a

    try:
        return converter[cls](a)
    except KeyError:
        for superclass in getmro(cls):
            try:
                return converter[superclass](a)
            except KeyError:
                continue

    try:
        return a._sympy_()
    except AttributeError:
        pass

    if not isinstance(a, basestring):
        for coerce in (float, int):
            try:
                return sympify(coerce(a))
            except (TypeError, ValueError, AttributeError, SympifyError):
                continue

    if strict:
        raise SympifyError(a)

    if iterable(a):
        try:
            return type(a)([sympify(x, locals=locals, convert_xor=convert_xor,
                rational=rational) for x in a])
        except TypeError:
            # Not all iterables are rebuildable with their type.
            pass
    if isinstance(a, dict):
        try:
            return type(a)([sympify(x, locals=locals, convert_xor=convert_xor,
                rational=rational) for x in a.iteritems()])
        except TypeError:
            # Not all iterables are rebuildable with their type.
            pass

    # At this point we were given an arbitrary expression
    # which does not inherit from Basic and doesn't implement
    # _sympy_ (which is a canonical and robust way to convert
    # anything to SymPy expression).
    #
    # As a last chance, we try to take "a"'s normal form via unicode()
    # and try to parse it. If it fails, then we have no luck and
    # return an exception
    try:
        a = unicode(a)
    except Exception, exc:
        raise SympifyError(a, exc)

    from sympy.parsing.sympy_parser import parse_expr, TokenError

    try:
        a = a.replace('\n', '')
        expr = parse_expr(a, locals or {}, rational, convert_xor)
    except (TokenError, SyntaxError):
        raise SympifyError('could not parse %r' % a)

    return expr

def _sympify(a):
    """Short version of sympify for internal usage for __add__ and __eq__
       methods where it is ok to allow some things (like Python integers
       and floats) in the expression. This excludes things (like strings)
       that are unwise to allow into such an expression.

       >>> from sympy import Integer
       >>> Integer(1) == 1
       True

       >>> Integer(1) == '1'
       False

       >>> from sympy import Symbol
       >>> from sympy.abc import x
       >>> x + 1
       x + 1

       >>> x + '1'
       Traceback (most recent call last):
           ...
       TypeError: unsupported operand type(s) for +: 'Symbol' and 'str'

       see: sympify
    """
    return sympify(a, strict=True)
