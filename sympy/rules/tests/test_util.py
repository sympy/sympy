from sympy.rules.util import treeexec

def test_treeexec():
    tree = ([3, 3], [4, 1], 2)
    assert treeexec(tree, {list: min, tuple: max}) == 3

    add = lambda *args: sum(args)
    mul = lambda *args: reduce(lambda a, b: a*b, args, 1)
    assert treeexec(tree, {list: add, tuple: mul}) == 60

