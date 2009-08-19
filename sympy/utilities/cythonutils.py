"""Helper module for cooperation with Cython. """

has_cython = True

try:
    import cython
except ImportError:
    has_cython = False

if has_cython:
    def cythonized(specs):
        arg_types = {}

        for spec in specs.split(','):
            arg_types[spec] = cython.int

        return cython.locals(**arg_types)
else:
    def cythonized(specs):
        return lambda f: f

