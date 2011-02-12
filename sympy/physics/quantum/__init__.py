__all__ = []

# Import all the name from modules listed below when using
# 'from sympy.physics.quantum import *'

modules = ['anticommutator', 'applyops', 'commutator', 'constants', 'dagger', \
    'hilbert', 'innerproduct', 'kronecker', 'operator', 'represent', 'state', \
    'tensorproduct'
]

# Prevent doctest running the code
if not __name__.endswith('__init__'):
    for m in modules:
        # Load all wanted modules...
        mod = __import__(__name__ + '.' + m, globals(), locals(), ['*'])
        for k in mod.__all__:
            # ...then store their exported functions in locals...
            locals()[k] = getattr(mod, k)
        # ...and let Python know they were exported.
        __all__.extend(mod.__all__)

