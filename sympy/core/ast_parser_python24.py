import compiler
import parser
from compiler.transformer import Transformer
from compiler.ast import CallFunc, Name, Const
from compiler.pycodegen import ExpressionCodeGenerator
import re

from basic import Basic
from symbol import Symbol

_is_integer = re.compile(r'\A\d+(l|L)?\Z').match

class SymPyTransformer(Transformer):
    def __init__(self, local_dict, global_dict):
        Transformer.__init__(self)

        self.symbol_class = 'Symbol'

        self.local_dict = local_dict
        self.global_dict = global_dict

    def atom_number(self, nodelist):
        n = Transformer.atom_number(self, nodelist)
        number, lineno = nodelist[0][1:]
        if _is_integer(number):
            n = Const(long(number), lineno)
            return CallFunc(Name('Integer'), [n])
        if number.endswith('j'):
            n = Const(complex(number), lineno)
            return CallFunc(Name('sympify'), [n])
        n = Const(number, lineno)
        return CallFunc(Name('Real'), [n])

    def atom_name(self, nodelist):
        name, lineno = nodelist[0][1:]

        if self.local_dict.has_key(name):
            name_obj = self.local_dict[name]
            return Const(name_obj, lineno=lineno)
        elif self.global_dict.has_key(name):
            name_obj = self.global_dict[name]

            if isinstance(name_obj, (Basic, type)) or callable(name_obj):
                return Const(name_obj, lineno=lineno)
        elif name in ['True', 'False']:
            return Const(eval(name), lineno=lineno)

        symbol_obj = Symbol(name)
        self.local_dict[name] = symbol_obj

        return Const(symbol_obj, lineno=lineno)

    def lambdef(self, nodelist):
        #this is python stdlib symbol, not SymPy symbol:
        from sympy import stdlib_symbol
        if nodelist[2][0] == stdlib_symbol.varargslist:
            names, defaults, flags = self.com_arglist(nodelist[2][1:])
        else:
            names = defaults = ()
            flags = 0

        lineno = nodelist[1][2]
        code = self.com_node(nodelist[-1])

        assert not defaults,`defaults` # sympy.Lambda does not support optional arguments

        if len(names) == 0:
            argument = ['x']
        else:
            argument = names

        def convert(x):
            return CallFunc(Name('sympify'), [Const(x)])

        argument = [ convert(arg) for arg in argument ]

        return CallFunc(Name('Lambda'), argument + [code])


class SymPyParser:
    def __init__(self, local_dict={}): #Contents of local_dict change, but it has proper effect only in global scope
        global_dict = {}
        exec 'from sympy import *' in global_dict

        self.r_transformer = SymPyTransformer(local_dict, global_dict)
        self.local_dict = local_dict
        self.global_dict = global_dict

    def parse_expr(self, ws_expression):
        expression = ws_expression.strip() #in case of "   x"
        ast_tree = parser.expr(expression)
        ast_tree = self.r_transformer.transform(ast_tree)

        compiler.misc.set_filename('<sympify>', ast_tree)
        code = ExpressionCodeGenerator(ast_tree).getCode()

        parsed_expr = eval(code, self.local_dict, self.global_dict) #Changed order to prefer sympy objects to  user defined

        return parsed_expr
