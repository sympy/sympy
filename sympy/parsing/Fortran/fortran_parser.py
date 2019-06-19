from lfortran.asr import asr
from lfortran.semantic import kinds
from lfortran.semantic.ast_to_asr import ast_to_asr
from lfortran.ast import src_to_ast, print_tree
from lfortran.asr.pprint import pprint_asr
from sympy.parsing.sympy_parser import (parse_expr, standard_transformations, implicit_multiplication_application)
from ast import *
import astor



class ASR2PyVisitor(asr.ASTVisitor):

    def __init__(self):
        self.py_ast = parse('')
        fix_missing_locations(self.py_ast)

    def visit_TranslationUnit(self, node):
        self._module = self.py_ast
        self._func = None

        for s in node.global_scope.symbols:
            sym = node.global_scope.symbols[s]
            self.visit(sym)
        for item in node.items:
            self.visit(item)

    def visit_Assignment(self, node):
        #TODO: Arithmatic Assignment
        if isinstance(node.target, asr.Variable):
            target = node.target
            value = node.value
            if (value.type == asr.Variable):
                new_node = Assign(targets = [Name(id = target.name, ctx = Store())], value = Name(id = value.name, ctx = Store()))
            else:
                if isinstance(value.op, asr.Add):
                    bin_op = Add()
                elif isinstance(value.op, asr.Sub):
                    bin_op = Sub()
                elif isinstance(value.op, asr.Div):
                    bin_op = Div()
                elif isinstance(value.op, asr.Mul):
                    bin_op = Mult()

                new_node = Assign(targets = [Name(id = target.name, ctx = Store())], value = BinOp(left = Name(id = value.left.name, ctx = Load()), op = bin_op, right = Name(id = value.right.name, ctz = Load())))
        else:
            raise NotImplementedError("Array")
        self.py_ast.body.append(new_node)
        fix_missing_locations(self.py_ast)

    def visit_BinOp(self, node):
        #TODO: Integer Binary Operations
        op = node.op
        lhs = node.left
        rhs = node.right
        #if isinstance(node.type, asr.Integer):
            #if isinstance(op, asr.Mul):
            #    new_node = Expr(value = BinOp(left = Num(lhs.n),
            #     op =Mult(), right = Num(rhs.num)))
            #elif isinstance(op, asr.Div):
            #    new_node = Expr(value = BinOp(left = Num(lhs.n),
            #     op =Div(), right = Num(rhs.num)))
            #elif isinstance(op, asr.Add):
            #    new_node = Expr(value = BinOp(left = Num(lhs.n),
            #     op =Add(), right = Num(rhs.num)))
            #elif isinstance(op, asr.Sub):
            #    new_node = Expr(value = BinOp(left = Num(lhs.n),
            #     op =Sub(), right = Num(rhs.num)))
            #else:
            #    raise NotImplementedError("Pow")

        if isinstance(node.type, asr.Variable) or isinstance(node.type, asr.Integer):
            if isinstance(op, asr.Mul):
                new_node = Expr(value = BinOp(left = Name(id = lhs.name, ctx = Load()),
                 op =Mult(), right = Name(id = rhs.name, ctx = Load())))
            elif isinstance(op, asr.Div):
                new_node = Expr(value = BinOp(left = Name(id = lhs.name, ctx = Load()),
                 op =Div(), right = Name(id = rhs.name, ctx = Load())))
            elif isinstance(op, asr.Add):
                new_node = Expr(value = BinOp(left = Name(id = lhs.name, ctx = Load()),
                 op =Add(), right = Name(id = rhs.name, ctx = Load())))
            elif isinstance(op, asr.Sub):
                new_node = Expr(value = BinOp(left = Name(id = lhs.name, ctx = Load()),
                 op =Sub(), right = Name(id = rhs.name, ctx = Load())))
            else:
                raise NotImplementedError("Pow")
        else:
            raise NotImplementedError()
        self.py_ast.body.append(new_node)
        fix_missing_locations(self.py_ast)

    def visit_Variable(self, node):
        new_node = Expr(value = Name(id = node.name, ctx = Load()))
        self.py_ast.body.append(new_node)
        fix_missing_locations(self.py_ast)

    def visit_sequence(self, seq):
        r = []
        if seq is not None:
            for node in seq:
                r.append(self.visit(node))
        self.py_ast.body.append(i for i in r)
        fix_missing_locations(self.py_ast)


    def visit_Num(self, node):
        py_ast.append(Num(n=node.n))
        fix_missing_locations(self.py_ast)

    def visit_Function(self, node):
        # TODO: Return statement, variable declaration
        fn_args =[]
        fn_name = node.name
        for arg_iter in node.args:
            fn_args.append(arg(arg = arg_iter.name, annotation=None))
        for i in node.body:
            fn_ast = call_visitor(i)
        fn_body = fn_ast.body
        #for sym in node.symtab.symbols:
        #    fn_body.append(Expr(value = Name(id = node.symtab.symbols[sym].name, ctx = Store())))
        fn_return = Name(id = node.return_var.name, ctx = Load())
        fn_body.append(Return(value = Name(id = fn_return, ctx = Store())))
        new_node = FunctionDef(name = fn_name, args = arguments(args = fn_args,vararg = None, kwarg =None, defaults = [], kw_defaults = []), decorator_list = [], body = fn_body)
        self.py_ast.body.append(new_node)
        fix_missing_locations(self.py_ast)

    def ret_ast(self):
        return self.py_ast

def call_visitor(fort_node):
    v = ASR2PyVisitor()
    v.visit(fort_node)
    res_ast = v.ret_ast()
    fix_missing_locations(res_ast)
    return res_ast


def src_to_sympy(src):
    a = ast_to_asr(src_to_ast(src,translation_unit=False))
    py_src = astor.to_source(call_visitor(a))
    print(py_src)
    #sym_src =  parse_expr("c = a+b",transformations=(standard_transformations + (implicit_multiplication_application,)))
    #print("SymPy Source Code")
    #print(sym_src)

t_src = """\
integer function f(a, b) result(r)
integer, intent(in) :: a, b
r = a + b
end function
"""
src_to_sympy(t_src)


def test_function():
        src1 = """\
        function abc(a ,b) result(inum2)
        integer :: i, j
        integer, intent(in) :: a,b
        integer, intent(out) :: inum2
        inum2 = a + b
        end function
        """
        res1 = """\
        def abc(a, b):
            inum = a + b
            return inum
        """
        py1 = src_to_sympy(src)
        assert res1.strip() == py1.strip()

def test_assignment():
    src2= """\
    integer :: a, b, c, d, e, f
    c = a + b
    d = a - b
    e = a * b
    f = a / b
    """
    res2 = """\
    a
    b
    c
    d
    e
    f
    c = a + b
    d = a - b
    e = a * b
    f = a / b
    """
    py2 = src_to_sympy(src2)
    assert res2.strip() == py2.strip()
