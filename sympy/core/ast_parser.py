"""
This module implements the functionality to take any Python expression as a
string and fix all numbers and other things before evaluating it,
thus

1/2

returns

Integer(1)/Integer(2)

We use the Python ast module for that, which is in python2.6 and later. It is
well documented at docs.python.org.

Some tips to understand how this works: use dump() to get a nice representation
of any node. Then write a string of what you want to get, e.g.
"Integer(1)", parse it, dump it and you'll see that you need to do
"Call(Name('Integer', Load()), [node], [], None, None)". You don't need to
bother with lineno and col_offset, just call fix_missing_locations() before
returning the node.

If the ast module is not available (python2.4 and 2.5), we use the old compiler
module.
"""

from sympy import Basic

try:
    import ast
    ast_enabled = True
except ImportError:
    ast_enabled = False

if ast_enabled:
    from ast import parse, NodeTransformer, Call, Name, Load, \
            fix_missing_locations, Str

    class Transform(NodeTransformer):

        def __init__(self, local_dict, global_dict):
            NodeTransformer.__init__(self)
            self.local_dict = local_dict
            self.global_dict = global_dict

        def visit_Num(self, node):
            if isinstance(node.n, int):
                return fix_missing_locations(Call(Name('Integer', Load()),
                        [node], [], None, None))
            elif isinstance(node.n, float):
                return fix_missing_locations(Call(Name('Real', Load()),
                    [node], [], None, None))
            return node

        def visit_Name(self, node):
            if node.id in self.local_dict:
                return node
            elif node.id in self.global_dict:
                name_obj = self.global_dict[node.id]

                if isinstance(name_obj, (Basic, type)) or callable(name_obj):
                    return node
            elif node.id in ['True', 'False']:
                return node
            return fix_missing_locations(Call(Name('Symbol', Load()),
                    [Str(node.id)], [], None, None))

        def visit_Lambda(self, node):
            if len(node.args.args) == 0:
                args = [Str("x")]
            else:
                args = node.args.args
            args = [self.visit(arg) for arg in args]
            body = self.visit(node.body)
            n = Call(Name('Lambda', Load()), args + [body], [], None, None)
            return fix_missing_locations(n)

def parse_expr(s, local_dict):
    """
    Converts the string "s" to a SymPy expression, in local_dict.

    It converts all numbers to Integers before feeding it to Python and
    automatically creates Symbols.
    """
    from sympify import SympifyError
    if ast_enabled:
        global_dict = {}
        exec 'from sympy import *' in global_dict
        try:
            a = parse(s.strip(), mode="eval")
        except SyntaxError:
            raise SympifyError("Cannot parse.")
        a = Transform(local_dict, global_dict).visit(a)
        e = compile(a, "<string>", "eval")
        return eval(e, local_dict, global_dict)
    else:
        # in python2.4 and 2.5, the "ast" module is not available, so we need
        # to use our old implementation:
        from ast_parser_python24 import SymPyParser
        try:
            return SymPyParser(local_dict=local_dict).parse_expr(s)
        except SyntaxError:
            raise SympifyError("sorry")
