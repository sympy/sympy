"""Module for SymPy containers

    (SymPy objects that store other SymPy objects)

    The containers implemented in this module are subclassed to Basic.
    They are supposed to work seamlessly within the SymPy framework.
"""
import sys

from sympy.core.basic import Basic
from sympy.core.sympify import sympify, converter
from sympy.utilities.iterables import iterable

class Tuple(Basic):
    """
    Wrapper around the builtin tuple object

    The Tuple is a subclass of Basic, so that it works well in the
    SymPy framework.  The wrapped tuple is available as self.args, but
    you can also access elements or slices with [:] syntax.

    >>> from sympy import symbols
    >>> from sympy.core.containers import Tuple
    >>> a, b, c, d = symbols('a b c d')
    >>> Tuple(a, b, c)[1:]
    (b, c)
    >>> Tuple(a, b, c).subs(a, d)
    (d, b, c)

    """

    def __new__(cls, *args, **assumptions):
        args = [ sympify(arg) for arg in args ]
        obj = Basic.__new__(cls, *args, **assumptions)
        return obj

    def __getitem__(self,i):
        if isinstance(i,slice):
            indices = i.indices(len(self))
            return Tuple(*[self.args[i] for i in range(*indices)])
        return self.args[i]

    def __len__(self):
        return len(self.args)

    def __contains__(self, item):
        return item in self.args

    def __iter__(self):
        return iter(self.args)

    def __add__(self, other):
        if isinstance(other, Tuple):
            return Tuple(*(self.args + other.args))
        elif isinstance(other, tuple):
            return Tuple(*(self.args + other))
        else:
            return NotImplemented

    def __radd__(self, other):
        if isinstance(other, Tuple):
            return Tuple(*(other.args + self.args))
        elif isinstance(other, tuple):
            return Tuple(*(other + self.args))
        else:
            return NotImplemented

    def __eq__(self, other):
        if isinstance(other, Basic):
            return super(Tuple, self).__eq__(other)
        return self.args == other

    def __ne__(self, other):
        if isinstance(other, Basic):
            return super(Tuple, self).__ne__(other)
        return self.args != other

    def __hash__(self):
        return hash(self.args)

    def _to_mpmath(self, prec):
        return tuple([a._to_mpmath(prec) for a in self.args])

    def __lt__(self, other):
        return self.args < other.args

    def __le__(self, other):
        return self.args <= other.args

converter[tuple] = lambda tup: Tuple(*tup)

def tuple_wrapper(method):
    """
    Decorator that converts any tuple in the function arguments into a Tuple.

    The motivation for this is to provide simple user interfaces.  The user can
    call a function with regular tuples in the argument, and the wrapper will
    convert them to Tuples before handing them to the function.

    >>> from sympy.core.containers import tuple_wrapper, Tuple
    >>> def f(*args):
    ...    return args
    >>> g = tuple_wrapper(f)

    The decorated function g sees only the Tuple argument:

    >>> g(0, (1, 2), 3)
    (0, (1, 2), 3)

    """
    def wrap_tuples(*args, **kw_args):
        newargs=[]
        for arg in args:
            if type(arg) is tuple:
                newargs.append(Tuple(*arg))
            else:
                newargs.append(arg)
        return method(*newargs, **kw_args)
    return wrap_tuples

class Dict(Basic):
    """
    Wrapper around the builtin dict object

    The Dict is a subclass of Basic, so that it works well in the
    SymPy framework.  Because it is immutable, it may be included
    in sets, but its values must all be given at instantiation and
    cannot be changed afterwards.  Otherwise it behaves identically
    to the Python dict.

    >>> from sympy import S
    >>> from sympy.core.containers import Dict

    >>> D = Dict({1: 'one', 2: 'two'})
    >>> for key in D:
    ...    if key == 1:
    ...        print key, D[key]
    1 one

    The args are sympified so the 1 and 2 are Integers and the values
    are Symbols. Queries automatically sympify args so the following work:

    >>> 1 in D
    True
    >>> D.has('one') # searches keys and values
    True
    >>> 'one' in D # not in the keys
    False
    >>> D[1]
    one

    """

    def __new__(cls, *args):
        if len(args)==1 and ((args[0].__class__ is dict) or
                             (args[0].__class__ is Dict)):
            items = [Tuple(k, v) for k, v in args[0].items()]
        elif iterable(args) and all(len(arg) == 2 for arg in args):
            items = [Tuple(k, v) for k, v in args]
        else:
            raise TypeError('Pass Dict args as Dict((k1, v1), ...) or Dict({k1: v1, ...})')
        elements = frozenset(items)
        obj = Basic.__new__(cls, elements)
        obj.elements = elements
        obj._dict = dict(items) # In case Tuple decides it wants to sympify
        return obj

    def __getitem__(self, key):
        """x.__getitem__(y) <==> x[y]"""
        return self._dict[sympify(key)]

    def __setitem__(self, key, value):
        raise NotImplementedError("SymPy Dicts are Immutable")

    @property
    def args(self):
        return tuple(self.elements)

    def items(self):
        '''D.items() -> list of D's (key, value) pairs, as 2-tuples'''
        return self._dict.items()

    def keys(self):
        '''D.keys() -> list of D's keys'''
        return self._dict.keys()

    def values(self):
        '''D.values() -> list of D's values'''
        return self._dict.values()

    def __iter__(self):
        '''x.__iter__() <==> iter(x)'''
        return iter(self._dict)

    def __len__(self):
        '''x.__len__() <==> len(x)'''
        return self._dict.__len__()

    def get(self, key, default=None):
        '''D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None.'''
        return self._dict.get(sympify(key), default)

    def __contains__(self, key):
        '''D.__contains__(k) -> True if D has a key k, else False'''
        return sympify(key) in self._dict

    def __lt__(self, other):
        return self.args < other.args

    @property
    def _sorted_args(self):
        from sympy.utilities import default_sort_key
        return sorted(self.args, key=default_sort_key)

if sys.version_info[:2] < (2, 6):
    from itertools import product
    ## {{{ http://code.activestate.com/recipes/576696/ (r5)
    import collections
    from weakref import proxy

    class Link(object):
        __slots__ = 'prev', 'next', 'key', '__weakref__'

    class OrderedSet(collections.MutableSet):
        'Set the remembers the order elements were added'
        # Big-O running times for all methods are the same as for regular sets.
        # The internal self.__map dictionary maps keys to links in a doubly linked list.
        # The circular doubly linked list starts and ends with a sentinel element.
        # The sentinel element never gets deleted (this simplifies the algorithm).
        # The prev/next links are weakref proxies (to prevent circular references).
        # Individual links are kept alive by the hard reference in self.__map.
        # Those hard references disappear when a key is deleted from an OrderedSet.

        def __init__(self, iterable=None):
            self.__root = root = Link()         # sentinel node for doubly linked list
            root.prev = root.next = root
            self.__map = {}                     # key --> link
            if iterable is not None:
                self |= iterable

        def __len__(self):
            return len(self.__map)

        def __contains__(self, key):
            return key in self.__map

        def add(self, key):
            # Store new key in a new link at the end of the linked list
            if key not in self.__map:
                self.__map[key] = link = Link()            
                root = self.__root
                last = root.prev
                link.prev, link.next, link.key = last, root, key
                last.next = root.prev = proxy(link)

        def discard(self, key):
            # Remove an existing item using self.__map to find the link which is
            # then removed by updating the links in the predecessor and successors.        
            if key in self.__map:        
                link = self.__map.pop(key)
                link.prev.next = link.next
                link.next.prev = link.prev

        def __iter__(self):
            # Traverse the linked list in order.
            root = self.__root
            curr = root.next
            while curr is not root:
                yield curr.key
                curr = curr.next

        def __reversed__(self):
            # Traverse the linked list in reverse order.
            root = self.__root
            curr = root.prev
            while curr is not root:
                yield curr.key
                curr = curr.prev

        def pop(self, last=True):
            if not self:
                raise KeyError('set is empty')
            key = next(reversed(self)) if last else next(iter(self))
            self.discard(key)
            return key

        def __repr__(self):
            if not self:
                return '%s()' % (self.__class__.__name__,)
            return '%s(%r)' % (self.__class__.__name__, list(self))

        def __eq__(self, other):
            if isinstance(other, OrderedSet):
                return len(self) == len(other) and list(self) == list(other)
            return not self.isdisjoint(other)

        def copy(self):
            return OrderedSet(list(self))
    ## end of http://code.activestate.com/recipes/576696/ }}}
else:
    ## {{{ http://code.activestate.com/recipes/576696/ (r5)
    import collections
    from weakref import proxy

    class Link(object):
        __slots__ = 'prev', 'next', 'key', '__weakref__'

    class OrderedSet(collections.MutableSet):
        'Set the remembers the order elements were added'
        # Big-O running times for all methods are the same as for regular sets.
        # The internal self.__map dictionary maps keys to links in a doubly linked list.
        # The circular doubly linked list starts and ends with a sentinel element.
        # The sentinel element never gets deleted (this simplifies the algorithm).
        # The prev/next links are weakref proxies (to prevent circular references).
        # Individual links are kept alive by the hard reference in self.__map.
        # Those hard references disappear when a key is deleted from an OrderedSet.

        def __init__(self, iterable=None):
            self.__root = root = Link()         # sentinel node for doubly linked list
            root.prev = root.next = root
            self.__map = {}                     # key --> link
            if iterable is not None:
                self |= iterable

        def __len__(self):
            return len(self.__map)

        def __contains__(self, key):
            return key in self.__map

        def add(self, key):
            # Store new key in a new link at the end of the linked list
            if key not in self.__map:
                self.__map[key] = link = Link()            
                root = self.__root
                last = root.prev
                link.prev, link.next, link.key = last, root, key
                last.next = root.prev = proxy(link)

        def discard(self, key):
            # Remove an existing item using self.__map to find the link which is
            # then removed by updating the links in the predecessor and successors.        
            if key in self.__map:        
                link = self.__map.pop(key)
                link.prev.next = link.next
                link.next.prev = link.prev

        def __iter__(self):
            # Traverse the linked list in order.
            root = self.__root
            curr = root.next
            while curr is not root:
                yield curr.key
                curr = curr.next

        def __reversed__(self):
            # Traverse the linked list in reverse order.
            root = self.__root
            curr = root.prev
            while curr is not root:
                yield curr.key
                curr = curr.prev

        def pop(self, last=True):
            if not self:
                raise KeyError('set is empty')
            key = next(reversed(self)) if last else next(iter(self))
            self.discard(key)
            return key

        def __repr__(self):
            if not self:
                return '%s()' % (self.__class__.__name__,)
            return '%s(%r)' % (self.__class__.__name__, list(self))

        def __eq__(self, other):
            if isinstance(other, OrderedSet):
                return len(self) == len(other) and list(self) == list(other)
            return not self.isdisjoint(other)

        def copy(self):
            return OrderedSet(list(self))
    ## end of http://code.activestate.com/recipes/576696/ }}}
