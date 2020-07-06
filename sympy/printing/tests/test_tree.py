from sympy.matrices.expressions import MatrixSymbol
from sympy.printing.tree import tree


def test_print_tree_MatAdd():
    A = MatrixSymbol('A', 3, 3)
    B = MatrixSymbol('B', 3, 3)
    # XXX After making the _assumptions cache deterministic,
    # this may not use dry run
    assert tree(A + B)


def test_print_tree_MatAdd_noassumptions():
    A = MatrixSymbol('A', 3, 3)
    B = MatrixSymbol('B', 3, 3)

    test_str = \
"""MatAdd: A + B
+-MatrixSymbol: A
| +-Str: A
| +-Integer: 3
| +-Integer: 3
+-MatrixSymbol: B
  +-Str: B
  +-Integer: 3
  +-Integer: 3
"""

    assert tree(A + B, assumptions=False) == test_str
