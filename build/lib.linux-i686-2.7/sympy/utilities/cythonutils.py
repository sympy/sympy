"""Helper module for cooperation with Cython. """

from sympy.external import import_module

cython = import_module('cython')

if cython:
    def cythonized(specs):
        arg_types = {}

        for spec in specs.split(','):
            arg_types[spec] = cython.int

        return cython.locals(**arg_types)
else:
    def cythonized(specs):
        return lambda f: f
