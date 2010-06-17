"""sympify -- convert objects SymPy internal format"""

from types import NoneType
from inspect import getmro

from core import BasicMeta

class SympifyError(ValueError):
    def __init__(self, expr, base_exc=None):
        self.expr = expr
        self.base_exc = base_exc
    def __str__(self):
        if self.base_exc is None:
            return "SympifyError: %r" % (self.expr,)

        return "Sympify of expression '%s' failed, because of exception being raised:\n%s: %s" % (self.expr, self.base_exc.__class__.__name__, str(self.base_exc))

sympy_classes = BasicMeta.all_classes

converter = {}

def sympify(a, locals=None, convert_xor=True, strict=False, rational=False):
    """
    Converts an arbitrary expression to a type that can be used inside sympy.

    For example, it will convert python ints into instance of sympy.Rational,
    floats into instances of sympy.Real, etc. It is also able to coerce symbolic
    expressions which inherit from Basic. This can be useful in cooperation
    with SAGE.

    It currently accepts as arguments:
       - any object defined in sympy (except matrices [TODO])
       - standard numeric python types: int, long, float, Decimal
       - strings (like "0.09" or "2e-19")
       - booleans, including `None` (will leave them unchanged)

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

    If the option `strict` is set to `True`, only the types for which an
    explicit conversion has been defined are converted. In the other
    cases, a SympifyError is raised.

    >>> sympify(True)
    True
    >>> sympify(True, strict=True)
    Traceback (most recent call last):
    ...
    SympifyError: SympifyError: True

    """
    try:
        cls = a.__class__
    except AttributeError:  #a is probably an old-style class object
        cls = type(a)
    if cls in sympy_classes:
        return a
    if cls in (bool, NoneType):
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

    if isinstance(a, (list, tuple, set)):
        return type(a)([sympify(x, locals=locals, convert_xor=convert_xor, rational=rational) for x in a])

    # At this point we were given an arbitrary expression
    # which does not inherit from Basic and doesn't implement
    # _sympy_ (which is a canonical and robust way to convert
    # anything to SymPy expression).
    #
    # As a last chance, we try to take "a"'s  normal form via unicode()
    # and try to parse it. If it fails, then we have no luck and
    # return an exception
    try:
        a = unicode(a)
    except Exception, exc:
        raise SympifyError(a, exc)

    # In the following,
    # o process the xor symbol ^ -> **
    #
    # o preprocess fractions because long int fractions aren't handled by ast_parser:
    #     >>> S('2222222222/11111111111')
    #     0L
    #     >>> S('2222222222*1/11111111111')
    #     2222222222/11111111111
    #     >>> S('222222222222/11111111111')
    #     20L
    #     >>> S('222222222222*1/11111111111')
    #     222222222222/11111111111
    # So if we replace / with *1/ it looks like things will work.
    # ...but we should really use a parser to make sure we aren't in a string
    # so the following aims to do basically this:
    # a = a.replace('/', '*1/').replace('*1/*1/', '//')
    import tokenize, re
    from tokenize import TokenError
    from StringIO import StringIO
    reps = []
    # we need a token holder and a position locator
    M = '\0'
    L = '\1'
    # check to be sure and come up with a better plan if necessary.
    assert M not in a and L not in a
    try:
        # clear quotes
        t = StringIO(a).readline
        for code, sym, start, stop, line in tokenize.generate_tokens(t):
            rep = sym
            if code == 3:
                a = a.replace(sym, M, 1)
                reps.append(rep)

        # clear repeating decimals, e.g. 3.4[31] -> (3+4*1/10+31*1/990)
        # Note: the *1/ is used to avoid the long integer problems that are
        # handled below.
        repeated = re.findall(r'(([^0-9.]*)([0-9]*)\.([0-9]*)\[([0-9]+)\])',a)
        for r in repeated:
            sym, pre, pre_d, post_d, repetend = r
            zeros = '0'*len(post_d)
            rep = '%s(%s+%s*1/1%s+%s*1/%s%s)' % (pre, pre_d or '0', post_d or '0', zeros, repetend, '9'*len(repetend), zeros)
            a = a.replace(sym, M, 1)
            reps.append(rep)

        # clear long int fractions, xor, and take this opportunity to
        # rationalize reals because if we don't do it now it will be
        # more expensive with nsimplify later
        t = StringIO(a).readline
        for code, sym, start, stop, line in tokenize.generate_tokens(t):
            rep = sym
            if rational and code == 2: # a number
                from sympy.core.numbers import Rational
                rep = "(%s)" % str(Rational(sym))
                rep = '*1/'.join(rep.split('/')) # in case we get a long int fraction
            elif code == 51: # an operator
                if sym == '/':
                    rep = '*1/'
                elif sym == '//':
                    pass
                elif sym == '^' and convert_xor:
                    rep = '**'
                else:
                    continue
            else:
                continue
            # if the line we are parsing is always 1 line long then we
            # could use the start info, but we'll just put a token in
            # and see how many Ms are to the left and then replace the
            # L with an M
            a = a.replace(sym, L, 1)
            reps.insert(a[:a.find(L)].count(M), rep)
            a = a.replace(L, M)

        a = a.split(M)
        anew = []
        for i, m in enumerate(reps):
            anew.append(a[i])
            anew.append(reps[i])
        anew.append(a[-1])
        a = ''.join(anew)
    except TokenError:
        raise SympifyError('%r had a tokenizing problem' % a)

    import ast_parser
    return ast_parser.parse_expr(a, locals or {})

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
    return sympify(a, strict=True)
