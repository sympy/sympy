from sympy.utilities.iterables import preorder_traversal

def pprint_nodes(subtrees):
    """
    Prettyprints systems of nodes.

    Example:

    >> print pprint_nodes(["a", "b1\nb2", "c"])
    +-a
    +-b1
    | b2
    +-c
    >>

    """
    def indent(s,type=1):
        x = s.split("\n")
        r = "+-%s\n"%x[0]
        for a in x[1:]:
            if a=="": continue
            if type==1:
                r += "| %s\n"%a
            else:
                r += "  %s\n"%a
        return r
    if len(subtrees)==0: return ""
    f="";
    for a in subtrees[:-1]:
        f += indent(a)
    f += indent(subtrees[-1],2)
    return f

def print_node(node):
    """
    Returns an information about the "node".

    This includes class name, string representation and assumptions.
    """
    s = "%s: %s\n" % (node.__class__.__name__, str(node))
    if len(node._assumptions) > 0:
        for a in node._assumptions:
            s += "%s: %s\n" % (a, node._assumptions[a])
    return s

def tree(node):
    """
    Returns a tree representation of "node" as a string.

    It uses print_node() together with pprint_nodes() on node.args recursively.

    See also: print_tree()
    """
    subtrees = []
    for arg in node.args:
        subtrees.append(tree(arg))
    s = print_node(node)+pprint_nodes(subtrees)
    return s

def print_tree(node):
    """
    Prints a tree representation of "node".

    In [1]: print_tree(x**2)
    Pow: x**2
    +-Symbol: x
    | comparable: False
    | noncommutative: False
    | commutative: True
    +-Integer: 2
      real: True
      comparable: True
      commutative: True
      infinitesimal: False
      nonzero: True
      unbounded: False
      noncommutative: False
      noninteger: False
      zero: False
      complex: True
      bounded: True
      rational: True
      integer: True
      imaginary: False
      finite: True
      irrational: False


    See also: tree()
    """
    print tree(node)
