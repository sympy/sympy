"""Tools for manipulation of expressions using paths. """

from sympy.core import Basic

class EPath(object):
    """
    Manipulate expressions using paths.

    EPath grammar in EBNF notation::

        literal   ::= /[A-Za-z_][A-Za-z_0-9]*/
        number    ::= /-?\d+/
        type      ::= literal
        attribute ::= literal "?"
        all       ::= "*"
        slice     ::= "[" number? (":" number? (":" number?)?)? "]"
        range     ::= all | slice
        query     ::= (type | attribute) ("|" (type | attribute))*
        selector  ::= range | query range?
        path      ::= "/" selector ("/" selector)*

    """

    __slots__ = ["_path", "_epath"]

    def __new__(cls, path):
        """Construct new EPath. """
        if isinstance(path, EPath):
            return path

        if not path:
            raise ValueError("empty EPath")

        _path = path

        if path[0] == '/':
            path = path[1:]
        else:
            raise NotImplementedError("non-root EPath")

        epath = []

        for selector in path.split('/'):
            selector = selector.strip()

            if not selector:
                raise ValueError("empty selector")

            index = 0

            for c in selector:
                if c.isalnum() or c == '_' or c == '|' or c == '?':
                    index += 1
                else:
                    break

            attrs = []
            types = []

            if index:
                elements = selector[:index]
                selector = selector[index:]

                for element in elements.split('|'):
                    element = element.strip()

                    if not element:
                        raise ValueError("empty element")

                    if element.endswith('?'):
                        attrs.append(element[:-1])
                    else:
                        types.append(element)

            span = None

            if selector == '*':
                pass
            else:
                if selector.startswith('['):
                    try:
                        i = selector.index(']')
                    except ValueError:
                        raise ValueError("expected ']', got EOL")

                    _span, span = selector[1:i], []

                    if ':' not in _span:
                        span = int(_span)
                    else:
                        for elt in _span.split(':', 3):
                            if not elt:
                                span.append(None)
                            else:
                                span.append(int(elt))

                        span = slice(*span)

                    selector = selector[i+1:]

                if selector:
                    raise ValueError("trailing characters in selector")

            epath.append((attrs, types, span))

        obj = object.__new__(cls)

        obj._path = _path
        obj._epath = epath

        return obj

    def __repr__(self):
        return "%s(%r)" % (self.__class__.__name__, self._path)

    def _get_ordered_args(self, expr):
        """Sort ``expr.args`` using printing order. """
        if expr.is_Add:
            return expr.as_ordered_terms()
        elif expr.is_Mul:
            return expr.as_ordered_factors()
        else:
            return expr.args

    def _hasattrs(self, expr, attrs):
        """Check if ``expr`` has any of ``attrs``. """
        for attr in attrs:
            if not hasattr(expr, attr):
                return False

        return True

    def _hastypes(self, expr, types):
        """Check if ``expr`` is any of ``types``. """
        _types = [ cls.__name__ for cls in expr.__class__.mro() ]
        return bool(set(_types).intersection(types))

    def _has(self, expr, attrs, types):
        """Apply ``_hasattrs`` and ``_hastypes`` to ``expr``. """
        if not (attrs or types):
            return True

        if attrs and self._hasattrs(expr, attrs):
            return True

        if types and self._hastypes(expr, types):
            return True

        return False

    def apply(self, func, expr):
        """
        Modify parts of an expression selected by a path.

        **Examples**

        >>> from sympy.simplify.epathtools import EPath
        >>> from sympy import sin, cos, E, symbols
        >>> x, y, z, t = symbols('x, y, z, t')

        >>> path = EPath("/*/[0]/Symbol")
        >>> expr = [((x, 1), 2), ((3, y), z)]
        >>> path.apply(lambda expr: expr**2, expr)
        [((x**2, 1), 2), ((3, y**2), z)]

        >>> path = EPath("/*/*/Symbol")
        >>> expr = t + sin(x + 1) + cos(x + y + E)
        >>> path.apply(lambda expr: 2*expr, expr)
        t + sin(2*x + 1) + cos(2*x + 2*y + E)

        """
        def _apply(path, func, expr):
            if not path:
                return func(expr)
            else:
                selector, path = path[0], path[1:]
                attrs, types, span = selector

                if isinstance(expr, Basic):
                    if not expr.is_Atom:
                        args, basic = self._get_ordered_args(expr), True
                    else:
                        return expr
                elif hasattr(expr, '__iter__'):
                    args, basic = expr, False
                else:
                    return expr

                args = list(args)

                if span is not None:
                    if type(span) == slice:
                        indices = range(*span.indices(len(args)))
                    else:
                        indices = [span]
                else:
                    indices = xrange(len(args))

                for i in indices:
                    try:
                        arg = args[i]
                    except IndexError:
                        continue

                    if self._has(arg, attrs, types):
                        args[i] = _apply(path, func, arg)

                if basic:
                    return expr.func(*args)
                else:
                    return expr.__class__(args)

        return _apply(self._epath, func, expr)

    def select(self, expr):
        """
        Retrieve parts of an expression selected by a path.

        **Examples**

        >>> from sympy.simplify.epathtools import EPath
        >>> from sympy import sin, cos, E
        >>> from sympy.abc import x, y, z, t

        >>> path = EPath("/*/[0]/Symbol")
        >>> expr = [((x, 1), 2), ((3, y), z)]

        >>> path.select(expr)
        [x, y]

        >>> path = EPath("/*/*/Symbol")
        >>> expr = t + sin(x + 1) + cos(x + y + E)

        >>> path.select(expr)
        [x, x, y]

        """
        result = []

        def _select(path, expr):
            if not path:
                result.append(expr)
            else:
                selector, path = path[0], path[1:]
                attrs, types, span = selector

                if isinstance(expr, Basic):
                    args = self._get_ordered_args(expr)
                elif hasattr(expr, '__iter__'):
                    args = expr
                else:
                    return

                if span is not None:
                    if type(span) == slice:
                        args = args[span]
                    else:
                        try:
                            args = [args[span]]
                        except IndexError:
                            return

                for arg in args:
                    if self._has(arg, attrs, types):
                        _select(path, arg)

        _select(self._epath, expr)
        return result

def epath(path, expr=None, func=None):
    """
    Manipulate parts of an expression selected by a path.

    This function allows to manipulate large nested expressions in single
    line of code, utilizing techniques to those applied in XML processing
    standards (e.g. XPath).

    If ``func`` is ``None``, :func:`epath` retrieves elements selected by
    the ``path``. Otherwise it applies ``func`` to each matching element.

    **Syntax**

    select all : "/*"
        Equivalent of ``for arg in args:``.
    select slice : "/[0]" | "/[1:5]" | "/[1:5:2]"
        Supports standard Python's slice syntax.
    select by type : "/list" | "/list|tuple"
        Emulates :func:`isinstance`.
    select by attribute : "/__iter__?"
        Emulates :func:`hasattr`.

    **Parameters**

    path : str | EPath
        A path as a string or a compiled EPath.
    expr : Basic | iterable
        An expression or a container of expressions.
    func : callable (optional)
        A callable that will be applied to matching parts.
    args : tuple (optional)
        Additional positional arguments to ``func``.
    kwargs : dict (optional)
        Additional keyword arguments to ``func``.

    **Examples**

    >>> from sympy.simplify.epathtools import epath
    >>> from sympy import sin, cos, E
    >>> from sympy.abc import x, y, z, t

    >>> path = "/*/[0]/Symbol"
    >>> expr = [((x, 1), 2), ((3, y), z)]

    >>> epath(path, expr)
    [x, y]
    >>> epath(path, expr, lambda expr: expr**2)
    [((x**2, 1), 2), ((3, y**2), z)]

    >>> path = "/*/*/Symbol"
    >>> expr = t + sin(x + 1) + cos(x + y + E)

    >>> epath(path, expr)
    [x, x, y]
    >>> epath(path, expr, lambda expr: 2*expr)
    t + sin(2*x + 1) + cos(2*x + 2*y + E)

    """
    _epath = EPath(path)

    if expr is None:
        return _epath
    if func is None:
        return _epath.select(expr)
    else:
        return _epath.apply(func, expr)

def eselect(expr, epath):
    """
    Select parts of an expression according to an EPath specification.

    **Examples**

    >>> from sympy.simplify.epathtools import eselect
    >>> from sympy import sin, cos, E, symbols
    >>> x, y, z, t = symbols("x, y, z, t")

    >>> expr = [((x, 1), 2), ((3, y), z)]
    >>> eselect(expr, "/*/[0]/Symbol")
    [x, y]

    >>> expr = t + sin(x + 1) + cos(x + y + E)
    >>> eselect(expr, "/*/*/Symbol")
    [x, x, y]

    """
    return EPath(epath).select(expr)

def eapply(func, expr, epath):
    """
    Modify an expression by applying a function to some of its parts.

    **Examples**

    >>> from sympy.simplify.epathtools import eapply
    >>> from sympy import sin, cos, E, symbols
    >>> x, y, z, t = symbols('x, y, z, t')

    >>> expr = [((x, 1), 2), ((3, y), z)]
    >>> eapply(lambda expr: expr**2, expr, "/*/[0]/Symbol")
    [((x**2, 1), 2), ((3, y**2), z)]

    >>> expr = t + sin(x + 1) + cos(x + y + E)
    >>> eapply(lambda expr: 2*expr, expr, "/*/*/Symbol")
    t + sin(2*x + 1) + cos(2*x + 2*y + E)

    """
    return EPath(epath).apply(func, expr)
