from sympy.printing.dot import (purestr, styleof, attrprint, dotnode,
        dotedges, dotprint)
from sympy import Symbol, Integer, Basic, Expr
from sympy.abc import x

def test_purestr():
    assert purestr(Symbol('x')) == "Symbol(x)"
    assert purestr(Basic(1, 2)) == "Basic(1, 2)"


def test_styleof():
    styles = [(Basic, {'color': 'blue', 'shape': 'ellipse'}),
              (Expr,  {'color': 'black'})]
    assert styleof(Basic(1), styles) == {'color': 'blue', 'shape': 'ellipse'}

    x = Symbol('x')
    assert styleof(x + 1, styles) == {'color': 'black', 'shape': 'ellipse'}

def test_attrprint():
    assert attrprint({'color': 'blue', 'shape': 'ellipse'}) == \
           '"color"="blue", "shape"="ellipse"'

def test_dotnode():

    assert dotnode(x) ==\
            '"Symbol(x)" ["color"="black", "label"="x", "shape"="ellipse"];'
    assert dotnode(x+2) == \
            '"Add(Integer(2), Symbol(x))" ["color"="black", "label"="Add", "shape"="ellipse"];'


def test_dotedges():
    assert sorted(dotedges(x+2)) == [
        '"Add(Integer(2), Symbol(x))" -> "Integer(2)";',
        '"Add(Integer(2), Symbol(x))" -> "Symbol(x)";'
        ]

def test_dotprint():
    text = dotprint(x+2)
    assert all(e in text for e in dotedges(x+2))
    assert all(n in text for n in map(dotnode, (x, Integer(2), x+2)))
    assert 'digraph' in text

def test_dotprint_depth():
    text = dotprint(3*x+2, depth=1)
    assert dotnode(3*x+2) in text
    assert dotnode(x) not in text

def test_Matrix_and_non_basics():
    from sympy import MatrixSymbol
    n = Symbol('n')
    assert dotprint(MatrixSymbol('X', n, n))
