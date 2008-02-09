"""sympify -- convert objects SymPy internal format"""
# from basic import Basic, BasicType, S
# from numbers  import Integer, Real
# from interval import Interval
import decimal

class SympifyError(ValueError):
    def __init__(self, expr, base_exc=None):
        self.expr = expr
        self.base_exc = base_exc
    def __str__(self):
        if self.base_exc is None:
            return "SympifyError: %s" % (self.expr,)

        return "Sympify of expression '%s' failed, because of exception being raised:\n%s: %s" % (self.expr, self.base_exc.__class__.__name__, str(self.base_exc))


def sympify(a, sympify_lists=False, locals= {}):
    """Converts an arbitrary expression to a type that can be used
       inside sympy. For example, it will convert python int's into
       instance of sympy.Rational, floats into intances of sympy.Real,
       etc. It is also able to coerce symbolic expressions which does
       inherit after Basic. This can be useful in cooperation with SAGE.

       It currently accepts as arguments:
           - any object defined in sympy (except maybe matrices [TODO])
           - standard numeric python types: int, long, float, Decimal
           - strings (like "0.09" or "2e-19")

       If sympify_lists is set to True then sympify will also accept
       lists, tuples and sets. It will return the same type but with
       all of the entries sympified.

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
    if isinstance(a, Basic):
        return a
    if isinstance(a, BasicType):
        return a
    elif isinstance(a, bool):
        raise NotImplementedError("bool support")
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
    elif (a.__class__ in [list,tuple]) and len(a) == 2:
        # isinstance causes problems in the issue #432, so we use .__class__
        return Interval(*a)
    elif isinstance(a, (list,tuple,set)) and sympify_lists:
        return type(a)([sympify(x, True) for x in a])
    elif hasattr(a, "_sympy_"):
        # the "a" implements _sympy_() method, that returns a SymPy
        # expression (by definition), so we just use it
        return a._sympy_()
    else:
        # XXX this is here because of cyclic-import issues
        from sympy.matrices import Matrix
        from sympy.polynomials import Polynomial

        if isinstance(a, Polynomial):
            return a
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

        try:
            import ast_parser
            return ast_parser.SymPyParser(local_dict=locals).parse_expr(a)
        except Exception, exc:
            raise SympifyError(a, exc)
    raise SympifyError("%r is NOT a valid SymPy expression" % a)

