"""
This module provides convenient functions to transform sympy expressions to
lambda functions which can be used to calculate numerical values very fast.
"""

from __future__ import division
from sympy.core.sympify import sympify

import inspect

# These are the namespaces the lambda functions will use.
MATH = {}
MPMATH = {}
NUMPY = {}
SYMPY = {}

# Mappings between sympy and other modules function names.
MATH_TRANSLATIONS = {
    "Abs":"fabs",
    "ceiling":"ceil",
    "E":"e",
    "ln":"log",
}

MPMATH_TRANSLATIONS = {
    "ceiling":"ceil",
    "chebyshevt":"chebyt",
    "chebyshevu":"chebyu",
    "E":"e",
    "I":"j",
    "ln":"log",
    #"lowergamma":"lower_gamma",
    "oo":"inf",
    #"uppergamma":"upper_gamma",
    "LambertW":"lambertw",
    "Matrix":"matrix",
    "conjugate":"conj",
}

NUMPY_TRANSLATIONS = {
    "acos":"arccos",
    "acosh":"arccosh",
    "arg":"angle",
    "asin":"arcsin",
    "asinh":"arcsinh",
    "atan":"arctan",
    "atan2":"arctan2",
    "atanh":"arctanh",
    "ceiling":"ceil",
    "E":"e",
    "im":"imag",
    "ln":"log",
    "Matrix":"matrix",
    "Max":"amax",
    "Min":"amin",
    "oo":"inf",
    "re":"real",
}

# Available modules:
MODULES = {
    "math":(MATH, MATH_TRANSLATIONS, ("from math import *",)),
    "mpmath":(MPMATH, MPMATH_TRANSLATIONS, ("from sympy.mpmath import *",)),
    "numpy":(NUMPY, NUMPY_TRANSLATIONS, ("from numpy import *",)),
    "sympy":(SYMPY, {}, ("from sympy.functions import *",
                         "from sympy.matrices import Matrix",
                         "from sympy import Integral, pi, oo, nan, zoo, E, I",
                         "from sympy.utilities.iterables import iff"))
}

def _import(module, reload="False"):
    """
    Creates a global translation dictionary for module.

    The argument module has to be one of the following strings: "math",
    "mpmath", "numpy", "sympy".
    These dictionaries map names of python functions to their equivalent in
    other modules.
    """
    if not module in MODULES:
        raise NameError("This module can't be used for lambdification.")
    namespace, translations, import_commands = MODULES[module]
    # Clear namespace or exit
    if namespace:
        # The namespace was already generated, don't do it again if not forced.
        if reload:
            namespace.clear()
        else:
            return

    # It's possible that numpy is not available.
    for import_command in import_commands:
        try:
            exec import_command in {}, namespace
        except ImportError:
            raise ImportError("Can't import %s with command %s" % (module, import_command))

    # Add translated names to namespace
    for sympyname, translation in translations.iteritems():
        namespace[sympyname] = namespace[translation]

def lambdify(args, expr, modules=None, printer=None, use_imps=True):
    """
    Returns a lambda function for fast calculation of numerical values.

    Usage:
    >>> from sympy import sqrt, sin
    >>> from sympy.utilities.lambdify import lambdify
    >>> from sympy.abc import x, y, z
    >>> f = lambdify(x, x**2)
    >>> f(2)
    4
    >>> f = lambdify((x,y,z), [z,y,x])
    >>> f(1,2,3)
    [3, 2, 1]
    >>> f = lambdify(x, sqrt(x))
    >>> f(4)
    2.0
    >>> f = lambdify((x,y), sin(x*y)**2)
    >>> f(0, 5)
    0.0

    If not specified differently by the user, Sympy functions are replaced as
    far as possible by either python-math, numpy (if available) or mpmath
    functions - exactly in this order.
    To change this behavior, the "modules" argument can be used.
    It accepts:
     - the strings "math", "mpmath", "numpy", "sympy"
     - any modules (e.g. math)
     - dictionaries that map names of sympy functions to arbitrary functions
     - lists that contain a mix of the arguments above. (Entries that are first
        in the list have higher priority)

    Examples:
    (1) Use one of the provided modules:
        >> f = lambdify(x, sin(x), "math")

        Attention: Functions that are not in the math module will throw a name
                   error when the lambda function is evaluated! So this would
                   be better:
        >> f = lambdify(x, sin(x)*gamma(x), ("math", "mpmath", "sympy"))

    (2) Use some other module:
        >> import numpy
        >> f = lambdify((x,y), tan(x*y), numpy)

        Attention: There are naming differences between numpy and sympy. So if
                   you simply take the numpy module, e.g. sympy.atan will not be
                   translated to numpy.arctan. Use the modified module instead
                   by passing the string "numpy":

        >> f = lambdify((x,y), tan(x*y), "numpy")
        >> f(1, 2)
        -2.18503986326
        >> from numpy import array
        >> f(array([1, 2, 3]), array([2, 3, 5]))
        [-2.18503986 -0.29100619 -0.8559934 ]

    (3) Use own dictionaries:
        >> def my_cool_function(x): ...
        >> dic = {"sin" : my_cool_function}
        >> f = lambdify(x, sin(x), dic)

        Now f would look like:
        >> lambda x: my_cool_function(x)

    Functions present in `expr` can also carry their own numerical
    implementations, in a callable attached to the ``_imp_``
    attribute.  Usually you attach this using the
    ``implemented_function`` factory:

    >>> from sympy.abc import x, y, z
    >>> from sympy.utilities.lambdify import lambdify, implemented_function
    >>> from sympy import Function
    >>> f = implemented_function(Function('f'), lambda x : x+1)
    >>> func = lambdify(x, f(x))
    >>> func(4)
    5

    ``lambdify`` always prefers ``_imp_`` implementations to
    implementations in other namespaces, unless the ``use_imps`` input
    parameter is False.
    """
    # If the user hasn't specified any modules, use what is available.
    if modules is None:
        # Use either numpy (if available) or python.math where possible.
        # XXX: This leads to different behaviour on different systems and
        #      might be the reason for irreproducible errors.
        try:
            _import("numpy")
            modules = ("math", "numpy", "mpmath", "sympy")
        except ImportError:
            modules = ("math", "mpmath", "sympy")

    # Get the needed namespaces.
    namespaces = []
    # First find any function implementations
    if use_imps:
        namespaces.append(_imp_namespace(expr))
    # Check for dict before iterating
    if isinstance(modules, dict) or not hasattr(modules, '__iter__'):
        namespaces.append(modules)
    else:
        namespaces += list(modules)
    # fill namespace with first having highest priority
    namespace = {}
    for m in namespaces[::-1]:
        buf = _get_namespace(m)
        namespace.update(buf)

    if hasattr(expr, "atoms") :
        #Try if you can extract symbols from the expression.
        #Move on if expr.atoms in not implemented.
        syms = expr.atoms()
        for term in syms:
            namespace.update({str(term): term})

    # Create lambda function.
    lstr = lambdastr(args, expr, printer=printer)

    return eval(lstr, namespace)

def _get_namespace(m):
    """
    This is used by _lambdify to parse its arguments.
    """
    if isinstance(m, str):
        _import(m)
        return MODULES[m][0]
    elif isinstance(m, dict):
        return m
    elif hasattr(m, "__dict__"):
        return m.__dict__
    else:
        raise TypeError("Argument must be either a string, dict or module but it is: %s" % m)

def lambdastr(args, expr, printer=None):
    """
    Returns a string that can be evaluated to a lambda function.

    >>> from sympy.abc import x, y, z
    >>> from sympy.utilities.lambdify import lambdastr
    >>> lambdastr(x, x**2)
    'lambda x: (x**2)'
    >>> lambdastr((x,y,z), [z,y,x])
    'lambda x,y,z: ([z, y, x])'

    """
    if printer is not None:
        if inspect.isfunction(printer):
            lambdarepr = printer
        else:
            if inspect.isclass(printer):
                lambdarepr = lambda expr: printer().doprint(expr)
            else:
                lambdarepr = lambda expr: printer.doprint(expr)
    else:
        #XXX: This has to be done here because of circular imports
        from sympy.printing.lambdarepr import lambdarepr

    # Transform everything to strings.
    expr = lambdarepr(expr)
    if isinstance(args, str):
        pass
    elif hasattr(args, "__iter__"):
        args = ",".join(str(a) for a in args)
    else:
        args = str(args)

    return "lambda %s: (%s)" % (args, expr)

def _imp_namespace(expr, namespace=None):
    """ Return namespace dict with function implementations

    We need to search for functions in anything that can be thrown at
    us - that is - anything that could be passed as `expr`.  Examples
    include sympy expressions, as well as tuples, lists and dicts that may
    contain sympy expressions.

    Parameters
    ----------
    expr : object
       Something passed to lambdify, that will generate valid code from
       ``str(expr)``.
    namespace : None or mapping
       Namespace to fill.  None results in new empty dict

    Returns
    -------
    namespace : dict
       dict with keys of implemented function names within `expr` and
       corresponding values being the numerical implementation of
       function

    Examples
    --------
    >>> from sympy.abc import x, y, z
    >>> from sympy.utilities.lambdify import implemented_function, _imp_namespace
    >>> from sympy import Function
    >>> f = implemented_function(Function('f'), lambda x : x+1)
    >>> g = implemented_function(Function('g'), lambda x : x*10)
    >>> namespace = _imp_namespace(f(g(x)))
    >>> sorted(namespace.keys())
    ['f', 'g']
    """
    # Delayed import to avoid circular imports
    from sympy.core.function import FunctionClass
    if namespace is None:
        namespace = {}
    # tuples, lists, dicts are valid expressions
    if isinstance(expr, (list, tuple)):
        for arg in expr:
            _imp_namespace(arg, namespace)
        return namespace
    elif isinstance(expr, dict):
        for key, val in expr.items():
            # functions can be in dictionary keys
            _imp_namespace(key, namespace)
            _imp_namespace(val, namespace)
        return namespace
    # sympy expressions may be Functions themselves
    func = getattr(expr, 'func', None)
    if isinstance(func, FunctionClass):
        imp = getattr(func, '_imp_', None)
        if not imp is None:
            name = expr.func.__name__
            if name in namespace and namespace[name] != imp:
                raise ValueError('We found more than one '
                                 'implementation with name '
                                 '"%s"' % name)
            namespace[name] = imp
    # and / or they may take Functions as arguments
    if hasattr(expr, 'args'):
        for arg in expr.args:
            _imp_namespace(arg, namespace)
    return namespace

def implemented_function(symfunc, implementation):
    """ Add numerical `implementation` to function `symfunc`

    `symfunc` can by a Function, or a name, in which case we make an
    anonymous function with this name.  The function is anonymous in the
    sense that the name is not unique in the sympy namespace.

    Parameters
    ----------
    symfunc : str or ``sympy.FunctionClass`` instance
       If str, then create new anonymous sympy function with this as
       name.  If `symfunc` is a sympy function, attach implementation to
       function
    implementation : callable
       numerical implementation of function for use in ``lambdify``

    Returns
    -------
    afunc : sympy.FunctionClass instance
       function with attached implementation

    Examples
    --------
    >>> from sympy.abc import x, y, z
    >>> from sympy.utilities.lambdify import lambdify, implemented_function
    >>> from sympy import Function
    >>> f = implemented_function(Function('f'), lambda x : x+1)
    >>> lam_f = lambdify(x, f(x))
    >>> lam_f(4)
    5
    """
    # Delayed import to avoid circular imports
    from sympy.core.function import UndefinedFunction
    # if name, create anonymous function to hold implementation
    if isinstance(symfunc, basestring):
        symfunc = UndefinedFunction(symfunc)
    # We need to attach as a method because symfunc will be a class
    symfunc._imp_ = staticmethod(implementation)
    return symfunc
