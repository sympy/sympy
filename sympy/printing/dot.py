from sympy import Basic, Expr, Symbol, Integer, Rational, Float

default_styles = [(Basic, {'color': 'blue', 'shape': 'ellipse'}),
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


def styleof(expr, styles=default_styles):
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

def attrprint(d, delimiter=', '):
    """ Print a dictionary of attributes

    >>> from sympy.printing.dot import attrprint
    >>> print attrprint({'color': 'blue', 'shape': 'ellipse'})
    "color"="blue", "shape"="ellipse"
    """
    return delimiter.join('"%s"="%s"'%item for item in sorted(d.items()))

def dotnode(expr, styles=default_styles, labelfunc=str):
    """ String defining a node

    >>> from sympy.printing.dot import dotnode
    >>> from sympy.abc import x
    >>> print dotnode(x)
    "Symbol(x)" ["color"="black", "label"="x", "shape"="ellipse"];
    """
    if hasattr(expr, 'dotnode'):
        return expr.dotnode(styles)
    style = styleof(expr, styles)

    if isinstance(expr, Basic) and not expr.is_Atom:
        label = str(expr.__class__.__name__)
    else:
        label = labelfunc(expr)
    style['label'] = label
    return '"%s" [%s];' % (purestr(expr), attrprint(style))


def dotedges(expr):
    """ List of strings for all expr->expr.arg pairs

    >>> from sympy.printing.dot import dotedges
    >>> from sympy.abc import x
    >>> for e in dotedges(x+2): print e
    "Add(Integer(2), Symbol(x))" -> "Integer(2)";
    "Add(Integer(2), Symbol(x))" -> "Symbol(x)";
    """
    if hasattr(expr, 'dotedges'):
        return expr.dotedges()
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

def dotprint(expr, styles=default_styles, atom=lambda x: not isinstance(x,
    Basic), maxdepth=None, labelfunc=str, **kwargs):
    """
    DOT description of a SymPy expression tree

    Options are

    ``styles``: Styles for different classes.  The default is ``[(Basic,
          {'color': 'blue', 'shape': 'ellipse'}), (Expr, {'color':
          'black'})]``

    ``atom``: Function used to determine if an arg is an atom.  The default is
          ``lambda x: not isinstance(x, Basic)``.  Another good choice is
          ``lambda x: not x.args``.

    ``maxdepth``: The maximum depth.  The default is None, meaning no limit.

    ``labelfunc``: How to label leaf nodes.  The default is ``str``.  Another
          good option is ``srepr``. For example with ``str``, the leaf nodes
          of ``x + 1`` are labeled, ``x`` and ``1``.  With ``srepr``, they
          are labeled ``Symbol('x')`` and ``Integer(1)``.

    Additional keyword arguments are included as styles for the graph.

    Examples
    ========

    >>> from sympy.printing.dot import dotprint
    >>> from sympy.abc import x
    >>> print dotprint(x+2)
    digraph{

    # Graph style
    "rankdir"="TD"

    #########
    # Nodes #
    #########

    "Symbol(x)" ["color"="black", "label"="x", "shape"="ellipse"];
    "Integer(2)" ["color"="black", "label"="2", "shape"="ellipse"];
    "Add(Integer(2), Symbol(x))" ["color"="black", "label"="Add", "shape"="ellipse"];

    #########
    # Edges #
    #########

    "Add(Integer(2), Symbol(x))" -> "Symbol(x)";
    "Add(Integer(2), Symbol(x))" -> "Integer(2)";
    }

    """

    graphstyle.update(kwargs)

    nodes = []
    edges = []
    def traverse(e, depth):
        nodes.append(dotnode(e, styles, labelfunc=labelfunc))
        if maxdepth and depth >= maxdepth:
            return
        edges.extend(dotedges(e))
        [traverse(arg, depth+1) for arg in e.args if not atom(arg)]
    traverse(expr, 0)

    return template%{'graphstyle': attrprint(graphstyle, delimiter='\n'),
                     'nodes': '\n'.join(sorted(set(nodes), key=len)),
                     'edges': '\n'.join(sorted(set(edges), key=len))}
