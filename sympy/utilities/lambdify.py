from __future__ import division

from sympy.core import sympify, vectorize

# These are the namespaces that lambda functions will use.
MATH = {}
MPMATH = {}
NUMPY = {}
SYMPY = {}

# Mappings between sympy and other modules function names.
MATH_TRANSLATIONS = {
    "ceiling":"ceil",
    "abs":"fabs",
    "ln":"log",
    "pi":"pi",
    "E":"e"
}

MPMATH_TRANSLATIONS = {
    "lowergamma":"lower_gamma",
    "uppergamma":"upper_gamma",
    "pi":"pi",
    "E":"e",
    "I":"j"
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
    "im":"imag",
    "ln":"log",
    "max_":"max",
    "min_":"min",
    "re":"real",
    "E":"e",
    "oo":"inf",
}

# Available modules:
MODULES = {
    "math":(MATH, MATH_TRANSLATIONS, "import math"),
    "mpmath":(MPMATH, MPMATH_TRANSLATIONS, "from sympy.thirdparty import mpmath"),
    "numpy":(NUMPY, NUMPY_TRANSLATIONS, "import numpy"),
    "sympy":(SYMPY, {}, "from sympy import functions")
}

def _import(modulename, reload="False"):
    """
    Creates the global translated dictionary that maps sympy function names to
    the corresponding functions.
    """
    if not MODULES.has_key(modulename):
        raise NameError, "This module can't be used for lambdification."
    namespace, translations, import_command = MODULES[modulename]
    # Clear namespace or exit
    if namespace:
        if reload:
            namespace.clear()
        else:
            # This namespace was already generated, don't do that again.
            return

    # It's possible that numpy is not available.
    try:
        exec import_command + " as module"
    except ImportError:
        raise ImportError, "Can't import %s with command %s"%(modulename,
                                                              import_command)

    # Add all names of the module to our working namspace
    namespace.update(vars(module))
    # Add translated names to namespace
    for sympyname, translation in translations.iteritems():
        namespace[sympyname] = getattr(module, translation)

def lambdify(args, expr, modulenames=None):
    """
    Returns a lambda function that can be used to calculate fast numerical
    values of expressions.

    Usage:

    >>> from sympy import symbols, sqrt, sin
    >>> x,y,z = symbols('xyz')
    >>> f = lambdify([x], x**2)
    >>> f(2)
    4
    >>> f = lambdify([x,y,z], [z,y,x])
    >>> f(1,2,3)
    [3, 2, 1]
    >>> f = lambdify([x], sqrt(x))
    >>> f(4)
    2.0
    >>> f = lambdify((x,y), sin(x*y)**2)
    >>> f(0, 5)
    0.0


    Sympy functions are replaced as far as possible by either numpy functions
    (if available) or python.math functions. For more precise values consider
    using lambdify_mpmath, which is slower but has arbitrary precission.
    The used modules can also be directly specified in "modulenames".
    Possible modulenames are:
      - math (=python.math)
      - mpmath
      - numpy
      - sympy
    Alternatively a dictionary can be passed instead of one modulename, that
    maps sympy function names to an arbitrary function.

    Example:

    >>> f = lambdify([x, y], sin(x), modulenames="math")
    >>> f(0, 5)
    0.0

    """
    # If the user specified the modules use those.
    if not modulenames is None:
        return _lambdify(args, expr, modulenames)

    return _lambdify(args, expr, ("math", "sympy"))

def _lambdify(args, expr, modulenames):
    """
    Returns a lambda function that takes args as arguments and uses functions
    defined in the modules.
    Possible modulenames are:
      - math (=python.math)
      - mpmath
      - numpy
      - sympy
    You can also pass more than one modulename - then functions that are not
    defined in one module, can be taken from another module.
    Alternatively a dictionary can be passed instead of one modulename, that
    maps sympy function names to an arbitrary function.
    The modules that are first in the list have higher priority.
    """
    if hasattr(modulenames, "__iter__"):
        modulenames = list(modulenames)
        modulenames.reverse()
        namespace = {}
        for m in modulenames:
            namespace.update(_get_namespace(m))
    else:
        namespace = _get_namespace(modulenames)

    lstr = lambdastr(args, expr)
    return eval(lstr, namespace)

def _get_namespace(m):
    """
    This is only used by _lambdify to parse it's arguments.
    """
    namespace = {}
    if isinstance(m, str):
        _import(m)
        namespace.update(MODULES[m][0])
    elif isinstance(m, dict):
        namespace.update(m)
    else:
        raise NotImplementedError
    return namespace

def lambdify_math(args, expr):
    """
    Can be used exactly like lambdify, but it only uses python.math and sympy
    functions.
    """
    return _lambdify(args, expr, ("math", "sympy"))

def lambdify_mpmath(args, expr):
    """
    Can be used exactly like lambdify, but it only uses mpmath and sympy
    functions.
    """
    return _lambdify(args, expr, ("mpmath", "sympy"))

def lambdify_numpy(args, expr, prefix="__convert__"):
    """
    Can be used exactly like lambdify, but it only uses numpy and sympy
    functions.
    """
    return _lambdify(args, expr, ("numpy", "sympy"))

def lambdastr(args, expr):
    """
    Returns a string that can be evaluated to a lambda function.
    lambdify() uses this.

    >>> from sympy import symbols
    >>> x,y,z = symbols('xyz')
    >>> lambdastr([x], x**2)
    'lambda x: (x**2)'
    >>> lambdastr([x,y,z], [z,y,x])
    'lambda x,y,z: ([z, y, x])'
    """
    # Transform everything to strings.
    expr = str(expr)
    if isinstance(args, str):
        pass
    elif hasattr(args, "__iter__"):
        args = ",".join(str(a) for a in args)
    elif args.is_Symbol:
        args = str(args)
    else:
        raise NotImplementedError

    # Set names of used variables back so that the expression can be used again.
    return "lambda %s: (%s)" % (args, expr)
