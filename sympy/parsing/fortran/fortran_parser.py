from sympy.codegen.ast import ( Variable, IntBaseType, FloatBaseType, String, Return, FunctionDefinition)
from sympy.external import import_module
from sympy.printing import pycode

lfortran = import_module(
    'lfortran',
    warn_not_installed=True
)

asr = lfortran.asr.asr
asr_mod = lfortran.asr
src_to_ast = lfortran.ast.src_to_ast
ast_to_asr = lfortran.semantic.ast_to_asr.ast_to_asr

"""
This module contains all the necessary Classes and Function used to Parse Fortran
code into SymPy expression

The module and its API are currently under development and experimental.
It is also dependent on LFortran for the ASR that is converted to SymPy syntax
which is also under development.
The module only supports the features currently supported by the LFortran ASR
which will be updated as the development of LFortran and this module progresses.

You might find unexpected bugs and exceptions while using the module, feel free
to report them to the SymPy Issue Tracker

The API for the module might also change while in development if better and more
effective ways are discovered for the process

Features Supported
==================

- Variable Declarations (integers and reals)
- Function Definitions


Notes
=====

The module currently depends on two external dependencies

LFortran : Required to parse Fortran source code into ASR


Refrences
=========

.. [1] https://github.com/sympy/sympy/issues
.. [2] https://gitlab.com/lfortran/lfortran
.. [3] https://docs.lfortran.org/
"""

class ASR2PyVisitor(asr.ASTVisitor):
    """
    Visitor Class for LFortran ASR

    It is a Visitor class derived from asr.ASRVisitor which visits all the nodes
    of the LFortran ASR and createds correspondiong Python AST node for each
    ASR node
    """
    def __init__(self):
        """Initialize the parser"""
        self.py_ast = []

    def visit_TranslationUnit(self, node):
        """
        Function to visit all the elements of the Translation Unit
        created by LFortran ASR
        """
        for s in node.global_scope.symbols:
            sym = node.global_scope.symbols[s]
            self.visit(sym)
        for item in node.items:
            self.visit(item)

    def visit_Assignment(self, node):
        """Visitor Function for Assignment

        Visits each Assignment is the LFortran ASR and creates corresponding
        assignment for SymPy.

        Notes
        =====

        The function currently only supports variable assignment and binary
        operation assignments of varying multitudes

        Raises
        ======

        NotImplementedError() when called for Numeric assignments or Arrays
        """
        #TODO: Arithmatic Assignment
        #if isinstance(node.target, asr.Variable):
        #    target = node.target
        #    value = node.value
        #    if isinstance(value, asr.Variable):
        #        new_node = Assign(
        #            targets = [
        #                Name(
        #                    id = target.name,
        #                    ctx = Store()
        #                )
        #            ],
        #            value = Name(
        #                id = value.name,
        #                ctx = Store()
        #            )
        #        )
        #    elif (type(value) == asr.BinOp):
        #        exp_ast = call_visitor_func(value)
        #        for expr in exp_ast.body:
        #            new_node = Assign(
        #                targets = [
        #                    Name(
        #                        id = target.name,
        #                        ctx = Store()
        #                    )
        #                ], value = expr
        #            )
        #    else:
        #        raise NotImplementedError("Numeric assignments not supported")
        #else:
        #    raise NotImplementedError("Arrays not supported")
        #self.py_ast.body.append(new_node)
        #fix_missing_locations(self.py_ast)
        raise NotImplementedError("Assignment not implemented in codegen AST")

    def visit_BinOp(self, node):
        """Visitor Function for Binary Operations

        Visits each binary operation present in the LFortran ASR like addition,
        substraction, multiplication,division and creates the operation node in
        the Python AST

        In case of more than one binary operations, the function calls the
        call_visitor_func() function on the child nodes pertaining to the
        binary operations recursively until all the operations have been included

        Notes
        =====

        The function currently only supports binary operations with Variables or
        other binary operation as nodes

        Raises
        ======

        NotImplementedError() when called for Numeric assignments
        """
        #TODO: Integer Binary Operations
        #op = node.op
        #lhs = node.left
        #rhs = node.right

        #if isinstance(op, asr.Add):
        #    bin_op = Add()
        #elif isinstance(op, asr.Sub):
        #    bin_op = Sub()
        #elif isinstance(op, asr.Div):
        #    bin_op = Div()
        #elif isinstance(op, asr.Mul):
        #    bin_op = Mult()

        #if (type(lhs) == asr.Variable):
        #    left_value = Name(
        #        id = lhs.name,
        #        ctx = Load()
        #    )
        #    if (type(rhs) == asr.Variable):
        #        new_node = BinOp(
        #            left = left_value,
        #            op = bin_op,
        #            right = Name(
        #                id = rhs.name,
        #                ctx = Load()
        #            )
        #        )
        #    elif (type(rhs) == asr.BinOp):
        #        r_exp_ast = call_visitor_func(rhs)
        #        for expr in r_exp_ast.body:
        #            new_node = BinOp(
        #                left = left_value,
        #                op = bin_op,
        #                right = expr
        #            )
        #    else:
        #        raise NotImplementedError("Numeric Assignments not supported")
        #else:
        #    l_exp_ast = call_visitor_func(lhs)
        #    for exp in l_exp_ast.body:
        #        left_value = exp
        #    if (type(rhs) == asr.Variable):
        #        new_node = BinOp(
        #            left = left_value,
        #            op = bin_op,
        #            right = Name(
        #                id = rhs.name,
        #                ctx = Load()
        #            )
        #        )
        #    elif(type(rhs) == asr.BinOp):
        #        r_exp_ast = call_visitor_func(rhs)
        #        for expr in r_exp_ast.body:
        #            new_node = BinOp(
        #                left = left_value,
        #                op = bin_op,
        #                right = expr
        #            )
        #    else:
        #        raise NotImplementedError("Numeric Assignments not supported")
        #self.py_ast.body.append(new_node)
        #fix_missing_locations(self.py_ast)
        raise NotImplementedError("Binary operations not implemented")

    def visit_Variable(self, node):
        """Visitor Function for Variable Declaration

        Visits each variable declaration present in the ASR and creates a
        Symbol declaration for each variable

        Notes
        =====
        The functions currently only supports the declaration of integer and
        real variables. Other data types are still under development.

        Raises
        ======

        NotImplementedError() when called for unsupported data types
        """
        if isinstance(node.type, asr.Integer):
            var_type = IntBaseType(String('integer'))
            value = 0
        elif isinstance(node.type, asr.Real):
            var_type = FloatBaseType(String('real'))
            value = 0.0
        else:
            raise NotImplementedError("Data type not supported")

        if not (node.intent == 'in'):
            new_node = Variable(
                node.name
            ).as_Declaration(
                type = var_type,
                value = value
            )               
            self.py_ast.append(new_node)

    def visit_Sequence(self, seq):
        """Visitor Function for code sequence

        Visits a code sequence/ block and calls the visitor function on all the
        children of the code block to create corresponding code in python
        """
        py_seq = []
        if seq is not None:
            for node in seq:
                expr = call_visitor(node)
                for elem in expr:
                    self.py_ast.append(elem)


    def visit_Num(self, node):
        """Visitor Function for Numbers in ASR

        This function is currently under development and will be updated
        with improvements in the LFortran ASR
        """
        #TODO:Numbers when the LFortran ASR is updated
        #py_ast.append(Num(n=node.n))
        #fix_missing_locations(self.py_ast)
        pass

    def visit_Function(self, node):
        """Visitor Function for function Definitions

        Visits each function definition present in the ASR and creates a
        function definition node in the Python AST with all the elements of the
        given function

        The functions declares all the variables required as SymPy symbols in
        the function before the function definition

        This function also the call_visior_function to parse the contents of
        the function body
        """
        # TODO: Return statement, variable declaration
        fn_args =[]
        fn_body = []
        fn_name = node.name
        for arg_iter in node.args:
            fn_args.append(
                Variable(
                    arg_iter.name
                )
            )
        for i in node.body:
            fn_ast = call_visitor(i)
        try:
            fn_body_expr = fn_ast
        except UnboundLocalError:
            fn_body_expr = []
        for sym in node.symtab.symbols:
            decl = call_visitor(node.symtab.symbols[sym])
            for symbols in decl:
                fn_body.append(symbols)
        for elem in fn_body_expr:
            fn_body.append(elem)
        fn_body.append(
            Return(
                Variable(
                    node.return_var.name
                )
            )
        )
        if isinstance(node.return_var.type, asr.Integer):
            ret_type = IntBaseType(String('integer'))
        elif isinstance(node.return_var.type, asr.Real):
            ret_type = FloatBaseType(String('real'))
        else:
            raise NotImplementedError("Data type not supported")
        new_node = FunctionDefinition(
                    return_type = ret_type,
                    name = fn_name,
                    parameters = fn_args,
                    body = fn_body
                )
        self.py_ast.append(new_node)

    def ret_ast(self):
        """Returns the AST"""
        return self.py_ast

def call_visitor(fort_node):
    """Calls the AST Visitor on the Module

    Thsi function is used to call the AST visitor for a program or module
    It imports all the required modules and calls the visit() function
    on the given node

    Parameters
    ==========

    fort_node : <lfortran.asr.asr.Assignment object> or <lfortran.asr.asr.BinOp object>
        The node for operation that the AST visitor is to be called upon

    Returns
    =======

    res_ast : list
        list of sympy AST Nodes
    """
    v = ASR2PyVisitor()
    v.visit(fort_node)
    res_ast = v.ret_ast()
    return res_ast

def src_to_sympy(src):
    """Wrapper function to convert the given Fortran source code to SymPy Expressions

    Parameters
    ==========

    src : string
        A string with the Fortran source code

    Returns
    =======

    py_src : string
        A string with the python source code compatible with SymPy
    """
    a_ast = src_to_ast(src,translation_unit=False)
    a = ast_to_asr(a_ast)
    py_src = call_visitor(a)
    return py_src
