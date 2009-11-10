"""
This module provides convenient functions to transform sympy expressions to
lambda functions which can be used to calculate numerical values very fast.
"""

from __future__ import division
from sympy.core.sympify import sympify

# These are the namespaces the lambda functions will use.
MATH = {}
MPMATH = {}
NUMPY = {}
SYMPY = {}

# Mappings between sympy and other modules function names.
MATH_TRANSLATIONS = {
    "abs":"fabs",
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
    "max_":"amax",
    "min_":"amin",
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
                         "from sympy import Integral"))
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

def lambdify(args, expr, modules=None):
    """
    Returns a lambda function for fast calculation of numerical values.

    Usage:
    >>> from sympy import sqrt, sin
    >>> from sympy.utilities import lambdify
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
    To change this behaviour, the "modules" argument can be used.
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

        Attention: There are naming diferences between numpy and sympy. So if
                   you simply take the numpy module, e.g. sympy.atan will not be
                   translated to numpy.arctan. Use the modified module instead
                   by passing the string "numpy".

    (3) Use own dictionaries:
        >> def my_cool_function(x): ...
        >> dic = {"sin" : my_cool_function}
        >> f = lambdify(x, sin(x), dic)

        Now f would look like:
        >> lambda x: my_cool_function(x)
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
    if isinstance(modules, dict): # Check for dict before "__iter__"
        namespace = _get_namespace(modules)
    elif hasattr(modules, "__iter__"):
        namespace = {}
        for m in modules:
            buf = _get_namespace(m)
            buf.update(namespace)
            namespace = buf
    else:
        namespace = _get_namespace(modules)

    # Create lambda function.
    lstr = lambdastr(args, expr)
    return eval(lstr, namespace)

def _get_namespace(m):
    """
    This is used by _lambdify to parse it's arguments.
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

def lambdastr(args, expr):
    """
    Returns a string that can be evaluated to a lambda function.

    >>> from sympy.abc import x, y, z
    >>> from sympy.utilities.lambdify import lambdastr
    >>> lambdastr(x, x**2)
    'lambda x: (x**2)'
    >>> lambdastr((x,y,z), [z,y,x])
    'lambda x,y,z: ([z, y, x])'

    """

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
