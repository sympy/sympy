from __future__ import unicode_literals, print_function

from sympy.codegen.ast import (
    Variable, IntBaseType, FloatBaseType, String,
    Integer, Float, FunctionPrototype, FunctionDefinition, FunctionCall,
    none, Return
    )
from sympy.printing import pycode
from sympy.external import import_module
import os
import sys

cin = import_module(
    'clang.cindex',
    warn_not_installed=True,
    __import__kwargs = {
        'fromlist':['cindex']
    }
)

"""
This module contains all the necessary Classes and Function used to Parse C and C++ code into SymPy expression
The module and its API are currently under development and experimental.
It is also dependent on Clang's AST adn Sympy's Cdegen AST that is converted to SymPy syntax
which might also bbe missing some features.
The module only supports the features currently supported by the Clang and codegen AST
which will be updated as the development of codegen AST and this module progresses.
You might find unexpected bugs and exceptions while using the module, feel free
to report them to the SymPy Issue Tracker
The API for the module might also change while in development if better and more
effective ways are discovered for the process

Features Supported
==================

- Variable Declarations (integers and reals)
- Assignment (using integer & floating literal and function calls)
- Function Definitions nad Declaration
- Function Calls
- Compound statements, Return statements

Notes
=====

The module is dependent on an external dependency which needs to be installed to use any features from the module.

Clang: The C and C++ compiler which is used to extract an AST from the provided C source code.

Refrences
=========

.. [1] https://github.com/sympy/sympy/issues
.. [2] https://clang.llvm.org/docs/
.. [3] https://clang.llvm.org/docs/IntroductionToTheClangAST.html
"""

class BaseParser(object):
    """Base Class for the C parser
    """
    def __init__(self):
        """Initializes the Base C parser creating a Clang AST index"""
        self.index = cin.Index.create()

    def diagnostics(self, out):
        """Diagostics function for the Clang AST"""
        for diag in self.tu.diagnostics:
            print('%s %s (line %s, col %s) %s' % (
                    {
                        4: 'FATAL',
                        3: 'ERROR',
                        2: 'WARNING',
                        1: 'NOTE',
                        0: 'IGNORED',
                    }[diag.severity],
                    diag.location.file,
                    diag.location.line,
                    diag.location.column,
                    diag.spelling
                ), file=out)

class CCodeConverter(BaseParser):
    """The Code Convereter for Clang AST


    """
    def __init__(self, name):
        """Initializes the code converter"""
        super(CCodeConverter, self).__init__()
        self._py_nodes = []

    def parse(self, filenames, flags):
        """Function to parse a file with C source code

        It takes the filename as an attribute and creates a Clang AST Translation Unit parsing the file.
        Then the transformation function is called on the transaltion unit, whose reults are collected into a list which is returned by the function.

        Parameters
        ==========

        filenames : string
            Path to the C file to be parsed

        flags: list
            Arguments to be passed to Clang while parsing the C code

        Returns
        =======

        py_nodes: list
            A list of sympy AST nodes

        """
        filename = os.path.abspath(filenames)
        self.tu = self.index.parse(
            filename,
            args=flags,
            options=cin.TranslationUnit.PARSE_DETAILED_PROCESSING_RECORD
        )
        for child in self.tu.cursor.get_children():
            if child.kind == cin.CursorKind.VAR_DECL:
                self._py_nodes.append(self.transform(child))
            elif (child.kind == cin.CursorKind.FUNCTION_DECL):
                self._py_nodes.append(self.transform(child))
            else:
                pass
        return self._py_nodes

    def transform(self, node):
        """Transformation Function for a Clang AST nodes

        It determines the kind of node and calss the respective transforation function for that node.

        Raises
        ======

        NotImplementedError : if the transformation for the provided node is not implemented

        """
        try:
            handler = getattr(self, 'transform_%s' % node.kind.name.lower())
        except AttributeError:
            print(
                "Ignoring node of type %s (%s)" % (
                    node.kind,
                    ' '.join(
                        t.spelling for t in node.get_tokens())
                    ),
                file=sys.stderr
            )
            handler = None
        if handler:
            result = handler(node)
            return result

    def transform_var_decl(self, node):
        """Transformation Function for Variable Declaration

        Used to create nodes for variable declarations and assignments with values or function call for the respective nodes in the clang AST

        Returns
        =======

        A variable node as Declaration, with the given value or 0 if the value is not provided

        Raises
        ======

        NotImplementedError : if called for data types not currently implemented

        Notes
        =====

        This function currently only supports basic Integer and Float data types

        """
        try:
            children = node.get_children()
            prev_child = None
            child = next(children)
            #ignoring namespace and type details for the variable
            while child.kind == cin.CursorKind.NAMESPACE_REF:
                prev_child = child
                child = next(children)

            while child.kind == cin.CursorKind.TYPE_REF:
                prev_child = child
                child = next(children)

            args = self.transform(child)
            # List in case of variable assignment, FunctionCall node in case of a funcion call
            if (child.kind == cin.CursorKind.INTEGER_LITERAL
                or child.kind == cin.CursorKind.UNEXPOSED_EXPR):
                return Variable(
                    node.spelling
                ).as_Declaration(
                    type = args[0],
                    value = args[1]
                )
            elif (child.kind == cin.CursorKind.CALL_EXPR):
                return Variable(
                    node.spelling
                ).as_Declaration(
                    value = args
                )
            else:
                raise NotImplementedError()

        except StopIteration:

            if (node.type.kind == cin.TypeKind.INT):
                type = IntBaseType(String('integer'))
                value = Integer(0)
            elif (node.type.kind == cin.TypeKind.FLOAT):
                type = FloatBaseType(String('real'))
                value = Float(0.0)
            else:
                raise NotImplementedError()
            return Variable(
                node.spelling
            ).as_Declaration(
                type = type,
                value = value
            )

    def transform_function_decl(self, node):
        """Transformation Function For Function Declaration

        Used to create nodes for function declarations and definitions for the respective nodes in the clang AST

        Returns
        =======

        function : Codegen AST node
            - FunctionPrototype node if function body is not present
            - FunctionDefinition node if the function body is present


        """
        token = node.get_tokens()
        c_ret_type = next(token).spelling
        if (c_ret_type == 'void'):
            ret_type = none
        elif(c_ret_type == 'int'):
            ret_type = type = IntBaseType(String('integer'))
        elif (c_ret_type == 'float'):
            ret_type = FloatBaseType(String('real'))
        else:
            raise NotImplementedError("Variable not yet supported")
        body = []
        param = []
        try:
            children = node.get_children()
            child = next(children)

            # If the node has any children, the first children will be the
            # return type and namespace for the function declaration. These
            # nodes can be ignored.
            while child.kind == cin.CursorKind.NAMESPACE_REF:
                child = next(children)

            while child.kind == cin.CursorKind.TYPE_REF:
                child = next(children)


            # Subsequent nodes will be the parameters for the function.
            try:

                while True:
                    decl = self.transform(child)
                    if (child.kind == cin.CursorKind.PARM_DECL):
                        param.append(decl)
                    elif (child.kind == cin.CursorKind.COMPOUND_STMT):
                        for val in decl:
                            body.append(val)
                    else:
                        body.append(decl)
                    child = next(children)
            except StopIteration:
                pass

        except StopIteration:
            pass

        if body == []:
            function = FunctionPrototype(
                return_type = ret_type,
                name = node.spelling,
                parameters = param
            )
        else:
            function = FunctionDefinition(
                return_type = ret_type,
                name = node.spelling,
                parameters = param,
                body = body
            )
        return function


    def transform_parm_decl(self, node):
        """Transformation function for Parameter Declaration

        Used to create parameter nodes for the required functions for the respective nodes in the clang AST

        Returns
        =======

        param : Codegen AST Node
            Variable node with the value nad type of the variable

        Raises
        ======

        ValueError if multiple children encountered in the parameter node

        """
        if (node.type.kind == cin.TypeKind.INT):
            type = IntBaseType(String('integer'))
            value = Integer(0)
        elif (node.type.kind == cin.TypeKind.FLOAT):
            type = FloatBaseType(String('real'))
            value = Float(0.0)
        try:
            children = node.get_children()
            child = next(children)

            # If there are any children, this will be a parameter
            # with a default value. The children will be the reference
            # to the default value.
            # If the default value is a non-primitive type, there will
            # be NAMESPACE_REF and TYPE_REF nodes; all but the last one
            # can be ignored.

            # Any namespace nodes can be stripped
            while child.kind in [cin.CursorKind.NAMESPACE_REF,
                                 cin.CursorKind.TYPE_REF,
                                 cin.CursorKind.TEMPLATE_REF]:
                child = next(children)

            # If there is a child, it is the default value of the parameter.
            args = self.transform(child)
            param = Variable(
                node.spelling
            ).as_Declaration(
                type = args[0],
                value = args[1]
            )
        except StopIteration:
            param = Variable(
                node.spelling
            ).as_Declaration(
                    type = type,
                    value = value
            )

        try:
            value = self.transform(next(children))
            raise ValueError("Can't handle multiple children on parameter")
        except StopIteration:
            pass

        return param

    def transform_integer_literal(self, node):
        """Transformation function for integer literal

        Used to get the value and type of the given integer literal.

        Returns
        =======

        val : list
            List with two arguments type and Value
            type contains the type of the integer
            value contains the value stored in the variable

        Notes
        =====

        Only Base Integer type supported for now

        """
        type = IntBaseType(String('integer'))
        try:
            value = next(node.get_tokens()).spelling
        except StopIteration:
            # No tokens
            value = Integer(node.literal)
        val = [type, value]
        return val

    def transform_floating_literal(self, node):
        """Transformation function for floating literal

        Used to get the value and type of the given floating literal.

        Returns
        =======

        val : list
            List with two arguments type and Value
            type contains the type of float
            value contains the value stored in the variable

        Notes
        =====

        Only Base Float type supported for now

        """
        type = FloatBaseType(String('real'))
        try:
            value = next(node.get_tokens()).spelling
        except (StopIteration, ValueError):
            # No tokens
            value = Float(node.literal)
        val = [type, value]
        return val


    def transform_string_literal(self, node):
        #TODO: No string type in AST
        #type =
        #try:
        #    value = next(node.get_tokens()).spelling
        #except (StopIteration, ValueError):
            # No tokens
        #    value = node.literal
        #val = [type, value]
        #return val
        pass

    def transform_character_literal(self, node):
        #TODO: No string Type in AST
        #type =
        #try:
        #    value = next(node.get_tokens()).spelling
        #except (StopIteration, ValueError):
            # No tokens
        #    value = node.literal
        #val = [type, value]
        #return val
        pass

    def transform_unexposed_decl(self,node):
        """Transformation function for unexposed declarations"""
        pass

    def transform_unexposed_expr(self, node):
        """Transformation function for unexposed expression

        Unexposed expressions are used to wrap float, double literals and expressions

        Returns
        =======

        expr : Codegen AST Node
            the result from the wrapped expression

        None : NoneType
            No childs are found for the node

        Raises
        ======

        ValueError if the expression contains multiple children

        """
        # Ignore unexposed nodes; pass whatever is the first
        # (and should be only) child unaltered.
        try:
            children = node.get_children()
            expr = self.transform(next(children))
        except StopIteration:
            return None

        try:
            next(children)
            raise ValueError("Unexposed expression has > 1 children.")
        except StopIteration:
            pass

        return expr

    def transform_decl_ref_expr(self, node):
        """Returns the name of the declaration reference"""
        return node.spelling

    def transform_call_expr(self, node):
        """Transformation function for a call expression

        Used to create function call nodes for the function calls present in the C code

        Returns
        =======

        FunctionCall : Codegen AST Node
            FunctionCall node with parameters if any parameters are present

        """
        param = []
        children = node.get_children()
        child = next(children)

        while child.kind == cin.CursorKind.NAMESPACE_REF:
            child = next(children)
        while child.kind == cin.CursorKind.TYPE_REF:
            child = next(children)

        first_child = self.transform(child)
        try:
            for child in children:
                arg = self.transform(child)
                if (child.kind == cin.CursorKind.INTEGER_LITERAL):
                    param.append(Integer(arg[1]))
                elif (child.kind == cin.CursorKind.FLOATING_LITERAL):
                    param.append(Float(arg[1]))
                else:
                    param.append(arg)
            return FunctionCall(first_child, param)

        except StopIteration:
            return FunctionCall(first_child)

    def transform_return_stmt(self, node):
        """Returns the Return Node for a return statement"""
        return Return(next(node.get_children()).spelling)

    def transform_compound_stmt(self, node):
        """Transformation function for compond statemets

        Returns
        =======

        expr : list
            list of Nodes for the expressions present in the statement

        None : NoneType
            if the compound statement is empty

        """
        try:
            expr = []
            children = node.get_children()
            for child in children:
                expr.append(self.transform(child))
        except StopIteration:
            return None
        return expr

    def transform_decl_stmt(self, node):
        """Transformation function for declaration statements

        These statements are used to wrap different kinds of declararions like variable or function declaration
        The function calls the transformer function for the child of the given node

        Returns
        =======

        statement : Codegen AST Node
            contains the node returned by the children node for the type of declaration

        Raises
        ======

        ValueError if multiple children present

        """
        try:
            children = node.get_children()
            statement = self.transform(next(children))
        except StopIteration:
            pass

        try:
            self.transform(next(children))
            raise ValueError("Don't know how to handle multiple statements")
        except StopIteration:
            pass

        return statement

def convert_c_file(filename):
    """Function for converting a C source file

    The function reads the source code present in the given file and parses it to give out sympy expressions

    Parameters
    ==========

    filename : string
        A string with the path to the C source file to be parsed

    Returns
    =======

    res_src : list
        List of Python expression strings

    Note
    ====

    The Parser and it's API are currently under development and the API and the supported features can be changed while in development.

    """
    res_src = []
    converter = CCodeConverter('output')
    p = converter.parse(filename, flags = [])
    for i in p:
        if not isinstance(i, FunctionPrototype):
            #Ignoring Function prototype
            res_src.append(pycode(i))
        if isinstance(i, FunctionDefinition):
            res_src.append(pycode(i))
    return res_src
