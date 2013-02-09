from sympy import Basic
from functools import partial

new = Basic.__new__

def is_leaf(x):
    return not isinstance(x, Basic) or x.is_Atom

def treeexec(tree, join, leaf=lambda x: x):
    """ Apply functions onto recursive containers (tree)

    join - a dictionary mapping container types to functions
      e.g. ``{list: chain, tuple: minimize}``

    Keys are containers/iterables.  Values are functions [a] -> a.

    Examples
    --------

    >>> from sympy.rules import treeexec
    >>> tree = ([3, 3], [4, 1])
    >>> treeexec(tree, {list: min, tuple: max})
    3

    >>> add = lambda *args: sum(args)
    >>> mul = lambda *args: reduce(lambda a, b: a*b, args, 1)
    >>> treeexec(tree, {list: add, tuple: mul})
    30
    """
    if type(tree) in join:
        return join[type(tree)](*map(partial(treeexec, join=join, leaf=leaf),
                                     tree))
    else:
        return leaf(tree)

