"""Useful utility decorators. """
from __future__ import annotations

from typing import TypeVar
import sys
import types
import inspect
from functools import wraps
from itertools import groupby

# Keep this import for backwards compatibility:
from sympy.external.mpmath import conserve_mpmath_dps # noqa: F401

from sympy.utilities.exceptions import sympy_deprecation_warning


T = TypeVar('T')
"""A generic type"""


def _rebuild_opaque(expr, args):
    """Rebuild an expression only if it retains its arguments' relationals."""
    from sympy.core.basic import Basic
    from sympy.core.relational import Relational

    args = tuple(args)
    result = expr.func(*args)
    if isinstance(expr, Basic) and expr.has(Relational):
        protected = set().union(*(arg.atoms(Relational)
            for arg in args if isinstance(arg, Basic)))
        if protected and (not isinstance(result, Basic) or not protected <= result.atoms(Relational)):
            return expr
    return result


def _relational_opaque(func=None, *, allow_relational=False):
    """Keep relational operands opaque to an algebraic transformation.

    This is for transformations that would otherwise use structural APIs
    such as ``rewrite``, ``replace`` or ``xreplace`` to alter relational
    operands. It preserves relations while transforming surrounding arguments.
    With ``allow_relational=True``, the function receives a whole relation
    so that it can dispatch to its existing relational implementation.
    In that mode, an existing ``_eval_<function name>`` hook also receives
    the whole object.
    Consecutive arithmetic arguments without relations are transformed
    together, preserving noncommutative order around protected operands.
    """
    if func is None:
        return lambda func: _relational_opaque(func, allow_relational=allow_relational)
    signature = inspect.signature(func)
    first_arg = next(iter(signature.parameters))
    eval_method = '_eval_' + func.__name__

    @wraps(func)
    def wrapper(*args, **kwargs):
        expr = args[0] if args else kwargs.get(first_arg)
        if getattr(expr, 'is_Atom', False):
            return func(*args, **kwargs)

        # Basic imports this module, so core imports must be deferred.
        from sympy.core.basic import Basic
        from sympy.core.expr import Expr
        from sympy.core.relational import Relational
        from sympy.core.sympify import sympify

        if not args and first_arg not in kwargs:
            return func(**kwargs)
        expr = sympify(expr)
        has_relational = isinstance(expr, Basic) and expr.has(Relational)
        if has_relational:
            signature.bind(*args, **kwargs)
        if args:
            args = args[1:]
        else:
            kwargs.pop(first_arg)

        def apply(expr):
            if allow_relational and getattr(expr, eval_method, None) is not None:
                return func(expr, *args, **kwargs)
            if isinstance(expr, Relational):
                if allow_relational:
                    return func(expr, *args, **kwargs)
                return expr
            if isinstance(expr, Basic) and (
                    not isinstance(expr, Expr) or expr.has(Relational)):
                if expr.is_Add or expr.is_Mul:
                    newargs = []
                    for has_relational, group in groupby(
                            expr.args, lambda a: a.has(Relational)):
                        if has_relational:
                            newargs.extend(apply(a) for a in group)
                        else:
                            newargs.append(func(expr.func(*group), *args, **kwargs))
                    newargs = tuple(newargs)
                else:
                    newargs = tuple(apply(a) for a in expr.args)
                if newargs != expr.args:
                    return _rebuild_opaque(expr, newargs)
                return expr
            return func(expr, *args, **kwargs)

        if has_relational:
            return apply(expr)
        return func(expr, *args, **kwargs)

    return wrapper


def threaded_factory(func, use_add):
    """A factory for ``threaded`` decorators. """
    from sympy.core import sympify
    from sympy.matrices import MatrixBase
    from sympy.utilities.iterables import iterable

    signature = inspect.signature(func)

    def apply(expr, *args, **kwargs):
        if getattr(expr, 'is_Relational', False):
            signature.bind(expr, *args, **kwargs)
            return expr
        return func(expr, *args, **kwargs)

    @wraps(func)
    def threaded_func(expr, *args, **kwargs):
        if isinstance(expr, MatrixBase):
            return expr.applyfunc(lambda f: apply(f, *args, **kwargs))
        elif iterable(expr):
            try:
                return expr.__class__([apply(f, *args, **kwargs) for f in expr])
            except TypeError:
                return expr
        else:
            expr = sympify(expr)

            if use_add and expr.is_Add:
                return expr.__class__(*[apply(f, *args, **kwargs) for f in expr.args])
            else:
                return apply(expr, *args, **kwargs)

    return threaded_func


def threaded(func):
    """Apply ``func`` to sub--elements of an object, including :class:`~.Add`.

    This decorator is intended to make it uniformly possible to apply a
    function to all elements of composite objects, e.g. matrices, lists, tuples
    and other iterable containers, or just expressions.

    Relationals are left unchanged: applying an arbitrary function to
    each side need not preserve the truth value or definedness of a relation.

    This version of :func:`threaded` decorator allows threading over
    elements of :class:`~.Add` class. If this behavior is not desirable
    use :func:`xthreaded` decorator.

    Functions using this decorator must have the following signature::

      @threaded
      def function(expr, *args, **kwargs):

    """
    return threaded_factory(func, True)


def xthreaded(func):
    """Apply ``func`` to sub--elements of an object, excluding :class:`~.Add`.

    This decorator is intended to make it uniformly possible to apply a
    function to all elements of composite objects, e.g. matrices, lists, tuples
    and other iterable containers, or just expressions.

    Relationals are left unchanged: applying an arbitrary function to
    each side need not preserve the truth value or definedness of a relation.

    This version of :func:`threaded` decorator disallows threading over
    elements of :class:`~.Add` class. If this behavior is not desirable
    use :func:`threaded` decorator.

    Functions using this decorator must have the following signature::

      @xthreaded
      def function(expr, *args, **kwargs):

    """
    return threaded_factory(func, False)


class no_attrs_in_subclass:
    """Don't 'inherit' certain attributes from a base class

    >>> from sympy.utilities.decorator import no_attrs_in_subclass

    >>> class A(object):
    ...     x = 'test'

    >>> A.x = no_attrs_in_subclass(A, A.x)

    >>> class B(A):
    ...     pass

    >>> hasattr(A, 'x')
    True
    >>> hasattr(B, 'x')
    False

    """
    def __init__(self, cls, f):
        self.cls = cls
        self.f = f

    def __get__(self, instance, owner=None):
        if owner == self.cls:
            if hasattr(self.f, '__get__'):
                return self.f.__get__(instance, owner)
            return self.f
        raise AttributeError


def doctest_depends_on(exe=None, modules=None, disable_viewers=None,
                       python_version=None, ground_types=None):
    """
    Adds metadata about the dependencies which need to be met for doctesting
    the docstrings of the decorated objects.

    ``exe`` should be a list of executables

    ``modules`` should be a list of modules

    ``disable_viewers`` should be a list of viewers for :func:`~sympy.printing.preview.preview` to disable

    ``python_version`` should be the minimum Python version required, as a tuple
    (like ``(3, 0)``)
    """
    dependencies = {}
    if exe is not None:
        dependencies['executables'] = exe
    if modules is not None:
        dependencies['modules'] = modules
    if disable_viewers is not None:
        dependencies['disable_viewers'] = disable_viewers
    if python_version is not None:
        dependencies['python_version'] = python_version
    if ground_types is not None:
        dependencies['ground_types'] = ground_types

    def skiptests():
        from sympy.testing.runtests import DependencyError, SymPyDocTests, PyTestReporter # lazy import
        r = PyTestReporter()
        t = SymPyDocTests(r, None)
        try:
            t._check_dependencies(**dependencies)
        except DependencyError:
            return True  # Skip doctests
        else:
            return False # Run doctests

    def depends_on_deco(fn):
        fn._doctest_depends_on = dependencies
        fn.__doctest_skip__ = skiptests

        if inspect.isclass(fn):
            fn._doctest_depdends_on = no_attrs_in_subclass(
                fn, fn._doctest_depends_on)
            fn.__doctest_skip__ = no_attrs_in_subclass(
                fn, fn.__doctest_skip__)
        return fn

    return depends_on_deco


def public(obj: T) -> T:
    """
    Append ``obj``'s name to global ``__all__`` variable (call site).

    By using this decorator on functions or classes you achieve the same goal
    as by filling ``__all__`` variables manually, you just do not have to repeat
    yourself (object's name). You also know if object is public at definition
    site, not at some random location (where ``__all__`` was set).

    Note that in multiple decorator setup (in almost all cases) ``@public``
    decorator must be applied before any other decorators, because it relies
    on the pointer to object's global namespace. If you apply other decorators
    first, ``@public`` may end up modifying the wrong namespace.

    Examples
    ========

    >>> from sympy.utilities.decorator import public

    >>> __all__ # noqa: F821
    Traceback (most recent call last):
    ...
    NameError: name '__all__' is not defined

    >>> @public
    ... def some_function():
    ...     pass

    >>> __all__ # noqa: F821
    ['some_function']

    """
    if isinstance(obj, types.FunctionType):
        ns = obj.__globals__
        name = obj.__name__
    elif isinstance(obj, (type(type), type)):
        ns = sys.modules[obj.__module__].__dict__
        name = obj.__name__
    else:
        raise TypeError("expected a function or a class, got %s" % obj)

    if "__all__" not in ns:
        ns["__all__"] = [name]
    else:
        ns["__all__"].append(name)

    return obj


def memoize_property(propfunc):
    """Property decorator that caches the value of potentially expensive
    ``propfunc`` after the first evaluation. The cached value is stored in
    the corresponding property name with an attached underscore."""
    attrname = '_' + propfunc.__name__
    sentinel = object()

    @wraps(propfunc)
    def accessor(self):
        val = getattr(self, attrname, sentinel)
        if val is sentinel:
            val = propfunc(self)
            setattr(self, attrname, val)
        return val

    return property(accessor)


def deprecated(message, *, deprecated_since_version,
               active_deprecations_target, stacklevel=3):
    '''
    Mark a function as deprecated.

    This decorator should be used if an entire function or class is
    deprecated. If only a certain functionality is deprecated, you should use
    :func:`~.warns_deprecated_sympy` directly. This decorator is just a
    convenience. There is no functional difference between using this
    decorator and calling ``warns_deprecated_sympy()`` at the top of the
    function.

    The decorator takes the same arguments as
    :func:`~.warns_deprecated_sympy`. See its
    documentation for details on what the keywords to this decorator do.

    See the :ref:`deprecation-policy` document for details on when and how
    things should be deprecated in SymPy.

    Examples
    ========

    >>> from sympy.utilities.decorator import deprecated
    >>> from sympy import simplify
    >>> @deprecated("""\
    ... The simplify_this(expr) function is deprecated. Use simplify(expr)
    ... instead.""", deprecated_since_version="1.1",
    ... active_deprecations_target='simplify-this-deprecation')
    ... def simplify_this(expr):
    ...     """
    ...     Simplify ``expr``.
    ...
    ...     .. deprecated:: 1.1
    ...
    ...        The ``simplify_this`` function is deprecated. Use :func:`simplify`
    ...        instead. See its documentation for more information. See
    ...        :ref:`simplify-this-deprecation` for details.
    ...
    ...     """
    ...     return simplify(expr)
    >>> from sympy.abc import x
    >>> simplify_this(x*(x + 1) - x**2) # doctest: +SKIP
    <stdin>:1: SymPyDeprecationWarning:
    <BLANKLINE>
    The simplify_this(expr) function is deprecated. Use simplify(expr)
    instead.
    <BLANKLINE>
    See https://docs.sympy.org/latest/explanation/active-deprecations.html#simplify-this-deprecation
    for details.
    <BLANKLINE>
    This has been deprecated since SymPy version 1.1. It
    will be removed in a future version of SymPy.
    <BLANKLINE>
      simplify_this(x)
    x

    See Also
    ========
    sympy.utilities.exceptions.SymPyDeprecationWarning
    sympy.utilities.exceptions.sympy_deprecation_warning
    sympy.utilities.exceptions.ignore_warnings
    sympy.testing.pytest.warns_deprecated_sympy

    '''
    decorator_kwargs = {"deprecated_since_version": deprecated_since_version,
               "active_deprecations_target": active_deprecations_target}
    def deprecated_decorator(wrapped):
        if hasattr(wrapped, '__mro__'):  # wrapped is actually a class
            class wrapper(wrapped):
                __doc__ = wrapped.__doc__
                __module__ = wrapped.__module__
                _sympy_deprecated_func = wrapped
                if '__new__' in wrapped.__dict__:
                    def __new__(cls, *args, **kwargs):
                        sympy_deprecation_warning(message, **decorator_kwargs, stacklevel=stacklevel)
                        return super().__new__(cls, *args, **kwargs)
                else:
                    def __init__(self, *args, **kwargs):
                        sympy_deprecation_warning(message, **decorator_kwargs, stacklevel=stacklevel)
                        super().__init__(*args, **kwargs)
            wrapper.__name__ = wrapped.__name__
        else:
            @wraps(wrapped)
            def wrapper(*args, **kwargs):
                sympy_deprecation_warning(message, **decorator_kwargs, stacklevel=stacklevel)
                return wrapped(*args, **kwargs)
            wrapper._sympy_deprecated_func = wrapped
        return wrapper
    return deprecated_decorator
