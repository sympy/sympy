"""Tools for managing evaluation contexts. """

from __future__ import print_function, division

from sympy.utilities.iterables import dict_merge
from sympy.polys.polyutils import PicklableWithSlots

__known_options__ = set(['frac', 'gens', 'wrt', 'sort', 'order', 'domain',
    'modulus', 'gaussian', 'extension', 'field', 'greedy', 'symmetric'])

__global_options__ = []

__template__ = """\
def %(option)s(_%(option)s):
    return Context(%(option)s=_%(option)s)
"""

for option in __known_options__:
    exec(__template__ % { 'option': option })


class Context(PicklableWithSlots):

    __slots__ = ['__options__']

    def __init__(self, dict=None, **options):
        if dict is not None:
            self.__options__ = dict_merge(dict, options)
        else:
            self.__options__ = options

    def __getattribute__(self, name):
        if name in __known_options__:
            try:
                return object.__getattribute__(self, '__options__')[name]
            except KeyError:
                return None
        else:
            return object.__getattribute__(self, name)

    def __str__(self):
        return 'Context(%s)' % ', '.join(
            [ '%s=%r' % (key, value) for key, value in self.__options__.items() ])

    def __and__(self, other):
        if isinstance(other, Context):
            return Context(**dict_merge(self.__options__, other.__options__))
        else:
            raise TypeError("a context manager expected, got %s" % other)

    def __enter__(self):
        raise NotImplementedError('global context')

    def __exit__(self, exc_type, exc_val, exc_tb):
        raise NotImplementedError('global context')


def register_context(func):
    def wrapper(self, *args, **kwargs):
        return func(*args, **dict_merge(self.__options__, kwargs))

    wrapper.__doc__ = func.__doc__
    wrapper.__name__ = func.__name__

    setattr(Context, func.__name__, wrapper)

    return func
