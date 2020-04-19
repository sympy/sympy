from sympy.external import import_module

lfortran = import_module('lfortran')

if lfortran:
    from sympy.codegen.ast import (Variable, IntBaseType, FloatBaseType, String,
                                   Return, FunctionDefinition, Assignment)
    from sympy.core import Add, Mul, Integer, Float
    from sympy import Symbol

    asr = lfortran.asr.asr

    """
    This module contains all the necessary Classes and Function used to Parse
    Fortran code into SymPy expression

    The module and its API are currently under development and experimental.
    It is also dependent on LFortran for the ASR that is converted to SymPy syntax
    which is also under development.
    The module only supports the features currently supported by the LFortran ASR
    which will be updated as the development of LFortran and this module progresses

    You might find unexpected bugs and exceptions while using the module, feel free
    to report them to the SymPy Issue Tracker

    The API for the module might also change while in development if better and
    more effective ways are discovered for the process

    Features Supported
    ==================

    - Variable Declarations (integers and reals)
    - Function Definitions
    - Assignments and Basic Binary Operations


    Notes
    =====

    The module depends on an external dependency

    LFortran : Required to parse Fortran source code into ASR


    Refrences
    =========

    .. [1] https://github.com/sympy/sympy/issues
    .. [2] https://gitlab.com/lfortran/lfortran
    .. [3] https://docs.lfortran.org/

    """


    class ASR2PyVisitor(asr.ASTVisitor):  # type: ignore
        """
        Visitor Class for LFortran ASR

        It is a Visitor class derived from asr.ASRVisitor which visits all the
        nodes of the LFortran ASR and creates corresponding AST node for each
        ASR node

        """

        def __init__(self):
            """Initialize the Parser"""
            self._py_ast = []

        def visit_TranslationUnit(self, node):
            """
            Function to visit all the elements of the Translation Unit
            created by LFortran ASR
            """

            for sym in node.global_scope.symbols.values():
                self._py_ast.append(self.visit(sym))
            for item in node.items:
                self._py_ast.append(self.visit(item))

        def visit_Assignment(self, node):
            """Visitor Function for Assignment

            Visits each Assignment is the LFortran ASR and creates corresponding
            assignment for SymPy.

            Notes
            =====

            The function currently only supports variable assignment and binary
            operation assignments of varying multitudes. Any type of numberS or
            array is not supported.

            Raises
            ======

            NotImplementedError() when called for Numeric assignments or Arrays

            """
            # TODO: Arithmatic Assignment
            if isinstance(node.target, asr.Variable):
                target = node.target
                value = node.value
                if isinstance(value, asr.Variable):
                    return Assignment(
                        Variable(
                                target.name
                            ),
                        Variable(
                                value.name
                            )
                    )
                elif (type(value) == asr.BinOp):
                    expr = self.visit(value)
                    return Assignment(
                        Variable(target.name),
                        expr
                    )
                else:
                    raise NotImplementedError("Numeric assignments not supported")
            else:
                raise NotImplementedError("Arrays not supported")

        def visit_BinOp(self, node):
            """Visitor Function for Binary Operations

            Visits each binary operation present in the LFortran ASR like addition,
            subtraction, multiplication, division and creates the corresponding
            operation node in SymPy's AST

            In case of more than one binary operations, the function calls the
            call_visitor() function on the child nodes of the binary operations
            recursively until all the operations have been processed.

            Notes
            =====

            The function currently only supports binary operations with Variables
            or other binary operations. Numerics are not supported as of yet.

            Raises
            ======

            NotImplementedError() when called for Numeric assignments

            """
            # TODO: Integer Binary Operations
            op = node.op
            lhs = node.left
            rhs = node.right

            if (type(lhs) == asr.Variable):
                left_value = Symbol(lhs.name)
            elif(type(lhs) == asr.BinOp):
                left_value = self.visit(lhs)
            else:
                raise NotImplementedError("Numbers Currently not supported")

            if (type(rhs) == asr.Variable):
                right_value = Symbol(rhs.name)
            elif(type(rhs) == asr.BinOp):
                right_value = self.visit(rhs)
            else:
                raise NotImplementedError("Numbers Currently not supported")

            if isinstance(op, asr.Add):
                return Add(left_value, right_value)
            elif isinstance(op, asr.Sub):
                return Add(left_value, -right_value)
            elif isinstance(op, asr.Div):
                return Mul(left_value, 1/right_value)
            elif isinstance(op, asr.Mul):
                return Mul(left_value, right_value)

        def visit_Variable(self, node):
            """Visitor Function for Variable Declaration

            Visits each variable declaration present in the ASR and creates a
            Symbol declaration for each variable

            Notes
            =====

            The functions currently only support declaration of integer and
            real variables. Other data types are still under development.

            Raises
            ======

            NotImplementedError() when called for unsupported data types

            """
            if isinstance(node.type, asr.Integer):
                var_type = IntBaseType(String('integer'))
                value = Integer(0)
            elif isinstance(node.type, asr.Real):
                var_type = FloatBaseType(String('real'))
                value = Float(0.0)
            else:
                raise NotImplementedError("Data type not supported")

            return Variable(
                node.name
            ).as_Declaration(
                type = var_type,
                value = value
            )

        def visit_Sequence(self, seq):
            """Visitor Function for code sequence

            Visits a code sequence/ block and calls the visitor function on all the
            children of the code block to create corresponding code in python

            """
            ls = []
            if seq is not None:
                for node in seq:
                    ls.append(self.visit(node))
            return ls

        def visit_Num(self, node):
            """Visitor Function for Numbers in ASR

            This function is currently under development and will be updated
            with improvements in the LFortran ASR

            """
            # TODO:Numbers when the LFortran ASR is updated
            # return Integer(node.n)
            pass

        def visit_Function(self, node):
            """Visitor Function for function Definitions

            Visits each function definition present in the ASR and creates a
            function definition node in the Python AST with all the elements of the
            given function

            The functions declare all the variables required as SymPy symbols in
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
            for sym in node.symtab.symbols.values():
                if sym.intent != 'in':
                    fn_body.append(self.visit(sym))
            for elem in node.body:
                fn_body.append(self.visit(elem))
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
            return FunctionDefinition(
                        return_type = ret_type,
                        name = fn_name,
                        parameters = fn_args,
                        body = fn_body
                    )

        def ret_ast(self):
            """Returns the AST nodes"""
            return self._py_ast
else:
    class ASR2PyVisitor():  # type: ignore
        def __init__(self, *args, **kwargs):
            raise ImportError('LFortran is not installed, cannot parse Fortran code')
