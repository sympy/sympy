"""
This module adds several functions for interactive source code inspection.
"""

import inspect

def source(object):
    """
    Prints the source code of a given object.
    """
    print inspect.getsource(object)
