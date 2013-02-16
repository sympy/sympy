import pydot
from sympy import Basic, Expr, Symbol, Integer, Rational, Float

styles = [(Basic, {'color': 'blue', 'shape': 'ellipse'}),
          (Expr,  {'color': 'black'})]


slotClasses = (Symbol, Integer, Rational, Float)
def purestr(x):
    """ A string that follows obj = type(obj)(*obj.args) exactly """
    if not isinstance(x, Basic):
        return str(x)
    if type(x) in slotClasses:
        args = [getattr(x, slot) for slot in x.__slots__]
    else:
        args = x.args
    return "%s(%s)"%(type(x).__name__, ', '.join(map(purestr, args)))


def styleof(expr, styles=styles):
    """ Merge style dictionaries in order

    >>> from sympy import Symbol, Basic, Expr
    >>> from sympy.printing.dot import styleof
    >>> styles = [(Basic, {'color': 'blue', 'shape': 'ellipse'}),
    ...           (Expr,  {'color': 'black'})]

    >>> styleof(Basic(1), styles)
    {'color': 'blue', 'shape': 'ellipse'}

    >>> x = Symbol('x')
    >>> styleof(x + 1, styles)  # this is an Expr
    {'color': 'black', 'shape': 'ellipse'}
    """
    style = dict()
    for typ, sty in styles:
        if isinstance(expr, typ):
            style.update(sty)
    return style


def dotnode(expr, styles=styles):
    """ String defining a node

    >>> from sympy.printing.dot import dotnode
    >>> from sympy.abc import x
    >>> dotnode(x)      # doctest: +SKIP
    '"Symbol(x)" ["color"="black", "shape"="ellipse", "label"="x"]'
    """
    style = styleof(expr, styles)

    if isinstance(expr, Basic) and not expr.is_Atom:
        label = str(expr.__class__.__name__)
    else:
        label = str(expr)
    style['label'] = label
    stylestring = ', '.join('"%s"="%s"'%item for item in style.items())
    return '"%s" [%s]' % (purestr(expr), stylestring)


def dotedges(expr):
    """ List of strings for all expr->expr.arg pairs

    >>> from sympy.printing.dot import dotedges
    >>> from sympy.abc import x
    >>> for e in dotedges(x+2): print e
    "Add(Integer(2), Symbol(x))" -> "Integer(2)";
    "Add(Integer(2), Symbol(x))" -> "Symbol(x)";
    """
    if not isinstance(expr, Basic):
        return []
    else:
        return ['"%s" -> "%s";'%(purestr(expr), purestr(arg)) for arg in expr.args]

template = \
"""digraph{

# Graph style
%(graphstyle)s

#########
# Nodes #
#########

%(nodes)s

#########
# Edges #
#########

%(edges)s
}"""

graphstyle = {'rankdir': 'TD'}

def dotprint(expr, **kwargs):
    """ DOT description of a SymPy expression tree

    >>> from sympy.printing.dot import dotprint
    >>> from sympy.abc import x
    >>> dotprint(x+2) # doctest: +SKIP
    """

    graphstyle.update(kwargs)

    stack = [expr]
    nodes = []
    edges = []
    while stack:
        e = stack.pop()
        nodes.append(dotnode(e))
        edges.extend(dotedges(e))
        if isinstance(e, Basic):
            stack.extend(e.args)

    graphstyles = ['%s=%s;'%item for item in graphstyle.items()]

    return template%{'graphstyle': '\n'.join(graphstyles),
                     'nodes': '\n'.join(nodes),
                     'edges': '\n'.join(edges)}
