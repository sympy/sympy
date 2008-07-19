"""Thirdparty Packages for internal use.
"""

import sys
import os

def import_thirdparty(lib):
    """
    Imports a thirdparty package "lib" by setting all paths correctly.

    At the moment, there is only the "pyglet" library, so we just put
    pyglet to sys.path temporarily, then import "lib" and then restore the path.
    With more packages, we'll just put them to sys.path as well.
    """

    seen = set()
    def new_import(name, globals={}, locals={}, fromlist=[]):
        if name in seen:
            return old_import(name, globals, locals, fromlist)
        seen.add(name)
        sys.path.insert(0, os.path.join(os.path.abspath(os.path.dirname( \
            __file__)), "pyglet"))
        try:
            m = old_import(name, globals, locals, fromlist)
        finally:
            del sys.path[0]
        return m
    import __builtin__
    old_import = __builtin__.__import__
    __builtin__.__import__ = new_import
    try:
        m = __import__(lib)
    finally:
        __builtin__.__import__ = old_import

    return m
