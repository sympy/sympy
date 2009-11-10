"""sympify -- convert objects SymPy internal format"""
# from basic import Basic, BasicType, S
# from numbers  import Integer, Real
import decimal

class SympifyError(ValueError):
    def __init__(self, expr, base_exc=None):
        self.expr = expr
        self.base_exc = base_exc
    def __str__(self):
        if self.base_exc is None:
            return "SympifyError: %s" % (self.expr,)

        return "Sympify of expression '%s' failed, because of exception being raised:\n%s: %s" % (self.expr, self.base_exc.__class__.__name__, str(self.base_exc))


def sympify(a, locals=None, convert_xor=True):
    """Converts an arbitrary expression to a type that can be used
       inside sympy. For example, it will convert python int's into
       instance of sympy.Rational, floats into instances of sympy.Real,
       etc. It is also able to coerce symbolic expressions which does
       inherit after Basic. This can be useful in cooperation with SAGE.

       It currently accepts as arguments:
           - any object defined in sympy (except maybe matrices [TODO])
           - standard numeric python types: int, long, float, Decimal
           - strings (like "0.09" or "2e-19")
           - booleans, including None (will leave them unchanged)

       If the argument is already a type that sympy understands, it will do
       nothing but return that value. This can be used at the begining of a
       function to ensure you are working with the correct type.

       >>> from sympy import *

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

    """
    # XXX instead of duplicating _sympify it would be better to call _sympify
    # directly from here, but a lot of SymPy still calls sympify (no '_') and
    # this will add unneccesary overhead.
    #
    # When everything settles, let's refactor this.
    #                                      -- kirr
    if locals is None:
        locals = {}
    if isinstance(a, (Basic, BasicType, bool)):
        return a
    if a is None:
        return a
    elif isinstance(a, (int, long)):
        return Integer(a)
    elif isinstance(a, (float, decimal.Decimal)):
        return Real(a)
    elif isinstance(a, complex):
        real, imag = map(sympify, (a.real, a.imag))
        return real + S.ImaginaryUnit * imag
    elif isinstance(a, (list,tuple,set)):
        return type(a)([sympify(x) for x in a])

    # let's see if 'a' implements conversion methods such as '_sympy_' or
    # '__int__', that returns a SymPy (by definition) or SymPy compatible
    # expression, so we just use it
    for methname, conv in [
            ('_sympy_',None),
            ('__float__', Real),
            ('__int__', Integer),
            ]:
        meth = getattr(a, methname, None)
        if meth is None:
            continue

        # we have to be careful -- calling Class.__int__() almost always is not
        # a good idea
        try:
            v = meth()
        except TypeError:
            continue

        if conv is not None:
            v = conv(v)

        return v

    else:
        # XXX this is here because of cyclic-import issues
        from sympy.matrices import Matrix

        if isinstance(a, Matrix):
            raise NotImplementedError('matrix support')

        if not isinstance(a, str):
            # At this point we were given an arbitrary expression
            # which does not inherit from Basic and doesn't implement
            # _sympy_ (which is a canonical and robust way to convert
            # anything to SymPy expression).
            #
            # As a last chance, we try to take "a"'s  normal form via str()
            # and try to parse it. If it fails, then we have no luck and
            # return an exception
            a = str(a)

        if convert_xor:
            a = a.replace('^','**')
        import ast_parser
        return ast_parser.parse_expr(a, locals)
    raise SympifyError("%r is NOT a valid SymPy expression" % a)


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
       1 + x

       >>> x + '1'
       Traceback (most recent call last):
           ...
       TypeError: unsupported operand type(s) for +: 'Symbol' and 'str'

       see: sympify
    """
    if isinstance(a, Basic):
        return a
    if isinstance(a, BasicType):
        return a
    elif isinstance(a, (int, long)):
        return Integer(a)
    elif isinstance(a, (float, decimal.Decimal)):
        return Real(a)
    elif isinstance(a, complex):
        real, imag = map(sympify, (a.real, a.imag))
        ireal, iimag = int(real), int(imag)

        if ireal + iimag*1j == a:
            return ireal + iimag*S.ImaginaryUnit
        return real + S.ImaginaryUnit * imag

    # let's see if 'a' implements conversion methods such as '_sympy_' or
    # '__int__', that returns a SymPy (by definition) or SymPy compatible
    # expression, so we just use it
    for methname, conv in [
            ('_sympy_',None),
            ('__float__', Real),
            ('__int__', Integer),
            ]:
        meth = getattr(a, methname, None)
        if meth is None:
            continue

        # we have to be careful -- calling Class.__int__() almost always is not
        # a good idea
        try:
            v = meth()
        except TypeError:
            continue

        if conv is not None:
            v = conv(v)

        return v

    raise SympifyError("%r is NOT a valid SymPy expression" % (a,))





from numbers import Integer, Real
from basic import Basic, BasicType, S
