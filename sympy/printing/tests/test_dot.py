from sympy.printing.dot import (purestr, get_purestr_short, styleof, attrprint, dotnode,
        dotedges, dotprint)
from sympy.core.basic import Basic
from sympy.core.expr import Expr
from sympy.core.numbers import (Float, Integer)
from sympy.core.singleton import S
from sympy.core.symbol import (Symbol, symbols)
from sympy.printing.repr import srepr
from sympy.abc import x


def test_purestr():
    assert purestr(Symbol('x')) == "Symbol('x')"
    assert purestr(Basic(S(1), S(2))) == "Basic(Integer(1), Integer(2))"
    assert purestr(Float(2)) == "Float('2.0', precision=53)"

    assert purestr(Symbol('x'), with_args=True) == ("Symbol('x')", ())
    assert purestr(Basic(S(1), S(2)), with_args=True) == \
            ('Basic(Integer(1), Integer(2))', ('Integer(1)', 'Integer(2)'))
    assert purestr(Float(2), with_args=True) == \
        ("Float('2.0', precision=53)", ())

def test_get_purestr_short():
    purestr_short = get_purestr_short()
    assert purestr_short(Symbol('x')) == "Symbol('x')"
    assert purestr_short(Basic(S(1), S(2))) == "node_0"
    assert purestr_short(Float(2)) == "Float('2.0', precision=53)"

    assert purestr_short(Symbol('x'), with_args=True) == ("Symbol('x')", ())
    assert purestr_short(Basic(S(1), S(2)), with_args=True) == ('node_0', ('Integer(1)', 'Integer(2)'))
    assert purestr_short(Float(2), with_args=True) == ("Float('2.0', precision=53)", ())


def test_styleof():
    styles = [(Basic, {'color': 'blue', 'shape': 'ellipse'}),
              (Expr,  {'color': 'black'})]
    assert styleof(Basic(S(1)), styles) == {'color': 'blue', 'shape': 'ellipse'}

    assert styleof(x + 1, styles) == {'color': 'black', 'shape': 'ellipse'}


def test_attrprint():
    assert attrprint({'color': 'blue', 'shape': 'ellipse'}) == \
           '"color"="blue", "shape"="ellipse"'

def test_dotnode():

    assert dotnode(x, repeat=False) == \
        '"Symbol(\'x\')" ["color"="black", "label"="x", "shape"="ellipse"];'
    assert dotnode(x+2, repeat=False) == \
        '"Add(Integer(2), Symbol(\'x\'))" ' \
        '["color"="black", "label"="Add", "shape"="ellipse"];', \
        dotnode(x+2,repeat=0)

    assert dotnode(x + x**2, repeat=False) == \
        '"Add(Symbol(\'x\'), Pow(Symbol(\'x\'), Integer(2)))" ' \
        '["color"="black", "label"="Add", "shape"="ellipse"];'
    assert dotnode(x + x**2, repeat=True) == \
        '"Add(Symbol(\'x\'), Pow(Symbol(\'x\'), Integer(2)))_()" ' \
        '["color"="black", "label"="Add", "shape"="ellipse"];'

def test_dotnode_shortnames():
    purestr_short = get_purestr_short()

    assert dotnode(x, repeat=False, nodenamefunc=purestr_short) == \
        '"Symbol(\'x\')" ["color"="black", "label"="x", "shape"="ellipse"];'
    assert dotnode(x+2, repeat=False, nodenamefunc=purestr_short) == \
        '"node_0" ["color"="black", "label"="Add", "shape"="ellipse"];', \
        dotnode(x+2,repeat=0, nodenamefunc=purestr_short)

    assert dotnode(x + x**2, repeat=False, nodenamefunc=purestr_short) == \
        '"node_1" ["color"="black", "label"="Add", "shape"="ellipse"];'
    assert dotnode(x + x**2, repeat=True, nodenamefunc=purestr_short) == \
        '"node_1_()" ["color"="black", "label"="Add", "shape"="ellipse"];'

def test_dotedges():
    assert sorted(dotedges(x+2, repeat=False)) == [
        '"Add(Integer(2), Symbol(\'x\'))" -> "Integer(2)";',
        '"Add(Integer(2), Symbol(\'x\'))" -> "Symbol(\'x\')";'
    ]
    assert sorted(dotedges(x + 2, repeat=True)) == [
        '"Add(Integer(2), Symbol(\'x\'))_()" -> "Integer(2)_(0,)";',
        '"Add(Integer(2), Symbol(\'x\'))_()" -> "Symbol(\'x\')_(1,)";'
    ]

def test_dotedges_shortnames():
    purestr_short = get_purestr_short()

    assert sorted(dotedges(x+2, repeat=False, nodenamefunc=purestr_short)) == [
        '"node_0" -> "Integer(2)";',
        '"node_0" -> "Symbol(\'x\')";'
    ]
    assert sorted(dotedges(x + 2, repeat=True, nodenamefunc=purestr_short)) == [
        '"node_0_()" -> "Integer(2)_(0,)";',
        '"node_0_()" -> "Symbol(\'x\')_(1,)";'
    ]

def test_dotprint():
    text = dotprint(x+2, repeat=False)
    assert text == \
"""digraph{

# Graph style
"ordering"="out"
"rankdir"="TD"

#########
# Nodes #
#########

"node_0" ["color"="black", "label"="Add", "shape"="ellipse"];
"Integer(2)" ["color"="black", "label"="2", "shape"="ellipse"];
"Symbol('x')" ["color"="black", "label"="x", "shape"="ellipse"];

#########
# Edges #
#########

"node_0" -> "Integer(2)";
"node_0" -> "Symbol('x')";
}"""

    text = dotprint(x+x**2, repeat=False)
    assert text == \
"""digraph{

# Graph style
"ordering"="out"
"rankdir"="TD"

#########
# Nodes #
#########

"node_0" ["color"="black", "label"="Add", "shape"="ellipse"];
"Symbol('x')" ["color"="black", "label"="x", "shape"="ellipse"];
"node_1" ["color"="black", "label"="Pow", "shape"="ellipse"];
"Symbol('x')" ["color"="black", "label"="x", "shape"="ellipse"];
"Integer(2)" ["color"="black", "label"="2", "shape"="ellipse"];

#########
# Edges #
#########

"node_0" -> "Symbol('x')";
"node_0" -> "node_1";
"node_1" -> "Symbol('x')";
"node_1" -> "Integer(2)";
}"""

    text = dotprint(x+x**2, repeat=True)
    assert text == \
"""digraph{

# Graph style
"ordering"="out"
"rankdir"="TD"

#########
# Nodes #
#########

"node_0_()" ["color"="black", "label"="Add", "shape"="ellipse"];
"Symbol('x')_(0,)" ["color"="black", "label"="x", "shape"="ellipse"];
"node_1_(1,)" ["color"="black", "label"="Pow", "shape"="ellipse"];
"Symbol('x')_(1, 0)" ["color"="black", "label"="x", "shape"="ellipse"];
"Integer(2)_(1, 1)" ["color"="black", "label"="2", "shape"="ellipse"];

#########
# Edges #
#########

"node_0_()" -> "Symbol('x')_(0,)";
"node_0_()" -> "node_1_(1,)";
"node_1_(1,)" -> "Symbol('x')_(1, 0)";
"node_1_(1,)" -> "Integer(2)_(1, 1)";
}"""

    text = dotprint(x**x, repeat=True)
    assert text == \
"""digraph{

# Graph style
"ordering"="out"
"rankdir"="TD"

#########
# Nodes #
#########

"node_0_()" ["color"="black", "label"="Pow", "shape"="ellipse"];
"Symbol('x')_(0,)" ["color"="black", "label"="x", "shape"="ellipse"];
"Symbol('x')_(1,)" ["color"="black", "label"="x", "shape"="ellipse"];

#########
# Edges #
#########

"node_0_()" -> "Symbol('x')_(0,)";
"node_0_()" -> "Symbol('x')_(1,)";
}"""

def test_dotprint_shortnames():
    purestr_short = get_purestr_short()

    text = dotprint(x+2, repeat=False, shortnodenames=True)
    assert all(e in text for e in dotedges(x+2, repeat=False, nodenamefunc=purestr_short))
    assert all(
        n in text for n in [dotnode(expr, repeat=False, nodenamefunc=purestr_short)
        for expr in (x, Integer(2), x+2)])
    assert 'digraph' in text

    purestr_short = get_purestr_short()
    text = dotprint(x+x**2, repeat=False, shortnodenames=True)
    assert all(e in text for e in dotedges(x+x**2, repeat=False, nodenamefunc=purestr_short))
    assert all(
        n in text for n in [dotnode(expr, repeat=False, nodenamefunc=purestr_short)
        for expr in (x, Integer(2), Integer(1)+Integer(1), x**2)]) #would also pass without int+int, because node_0 would then be x**2
    assert 'digraph' in text

    purestr_short = get_purestr_short()
    text = dotprint(x+x**2, repeat=True, shortnodenames=True)
    assert all(e in text for e in dotedges(x+x**2, repeat=True, nodenamefunc=purestr_short))
    assert all(
        n in text for n in [dotnode(expr, pos=(), nodenamefunc=purestr_short)
        for expr in [x + x**2]])

    purestr_short = get_purestr_short()
    text = dotprint(x**x, repeat=True, shortnodenames=True)
    assert all(e in text for e in dotedges(x**x, repeat=True, nodenamefunc=purestr_short))
    assert all(
        n in text for n in [dotnode(x, pos=(0,), nodenamefunc=purestr_short),
                           dotnode(x, pos=(1,), nodenamefunc=purestr_short)])
    assert 'digraph' in text

def test_dotprint_depth():
    text = dotprint(3*x+2, depth=1)
    assert dotnode(3*x+2) not in text
    assert dotnode(x) not in text
    text = dotprint(3*x+2)
    assert "depth" not in text

def test_dotprint_depth_shortnames():
    text = dotprint(3*x+2, depth=1, shortnodenames=True)
    assert "node_0" in text
    assert dotnode(x) not in text
    text = dotprint(3*x+2, shortnodenames=True)
    assert "depth" not in text

def test_Matrix_and_non_basics():
    from sympy.matrices.expressions.matexpr import MatrixSymbol
    n = Symbol('n')
    text = dotprint(MatrixSymbol('X', n, n))
    assert text == \
"""digraph{

# Graph style
"ordering"="out"
"rankdir"="TD"

#########
# Nodes #
#########

"node_0_()" ["color"="black", "label"="MatrixSymbol", "shape"="ellipse"];
"Str('X')_(0,)" ["color"="blue", "label"="X", "shape"="ellipse"];
"Symbol('n')_(1,)" ["color"="black", "label"="n", "shape"="ellipse"];
"Symbol('n')_(2,)" ["color"="black", "label"="n", "shape"="ellipse"];

#########
# Edges #
#########

"node_0_()" -> "Str('X')_(0,)";
"node_0_()" -> "Symbol('n')_(1,)";
"node_0_()" -> "Symbol('n')_(2,)";
}"""

def test_labelfunc():
    text = dotprint(x + 2, labelfunc=srepr)
    assert "Symbol('x')" in text
    assert "Integer(2)" in text


def test_commutative():
    x, y = symbols('x y', commutative=False)
    assert dotprint(x + y) == dotprint(y + x)
    assert dotprint(x*y) != dotprint(y*x)
