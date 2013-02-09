from sympy.rules.util import treeexec

def test_treeexec():
    tree = ([3, 3], [4, 1], 2)
    assert treeexec(tree, {list: min, tuple: max}) == 3

    add = lambda *args: sum(args)
    mul = lambda *args: reduce(lambda a, b: a*b, args, 1)
    assert treeexec(tree, {list: add, tuple: mul}) == 60

def test_treeexec_leaf():
    assert treeexec(3, {'leaf': lambda x: x**2}) == 9
    tree = ([3, 3], [4, 1], 2)
    treep1 = ([4, 4], [5, 2], 3)
    assert treeexec(tree, {list: min, tuple: max, 'leaf': lambda x: x+1}) == \
           treeexec(treep1, {list: min, tuple: max})

