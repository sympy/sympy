"""
This module adds several functions for interactive source code inspection.
"""

import inspect

def source(object):
    """
    Prints the source code of a given object.
    """
    print 'In file: %s' % inspect.getsourcefile(object)
    print inspect.getsource(object)

def get_class(lookup_view):
    """
    Convert a string version of a class name to the object.

    For example, get_class('sympy.core.Basic') will return
    class Basic located in module sympy.core
    """
    if isinstance(lookup_view, str):
        # Bail early for non-ASCII strings (they can't be functions).
        lookup_view = lookup_view.encode('ascii')
        mod_name, func_name = get_mod_func(lookup_view)
        if func_name != '':
            lookup_view = getattr(__import__(mod_name, {}, {}, ['']), func_name)
            if not callable(lookup_view):
                raise AttributeError("'%s.%s' is not a callable." % (mod_name, func_name))
    return lookup_view

def get_mod_func(callback):
    """
    splits the string path to a class into a string path to the module
    and the name of the class. For example:
        >>> get_mod_func('sympy.core.basic.Basic')
        ('sympy.core.basic', 'Basic')

    """
    try:
        dot = callback.rindex('.')
    except ValueError:
        return callback, ''
    return callback[:dot], callback[dot+1:]
