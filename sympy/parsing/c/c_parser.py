from __future__ import unicode_literals, print_function
from sympy.external import import_module
import os

cin = import_module('clang.cindex', import_kwargs = {'fromlist': ['cindex']})

"""
This module contains all the necessary Classes and Function used to Parse C and
C++ code into SymPy expression
The module serves as a backend for SymPyExpression to parse C code
It is also dependent on Clang's AST and Sympy's Codegen AST.
The module only supports the features currently supported by the Clang and
codegen AST which will be updated as the development of codegen AST and this
module progresses.
You might find unexpected bugs and exceptions while using the module, feel free
to report them to the SymPy Issue Tracker

Features Supported
==================

- Variable Declarations (integers and reals)
- Assignment (using integer & floating literal and function calls)
- Function Definitions nad Declaration
- Function Calls
- Compound statements, Return statements

Notes
=====

The module is dependent on an external dependency which needs to be installed
to use the features of this module.

Clang: The C and C++ compiler which is used to extract an AST from the provided
C source code.

Refrences
=========

.. [1] https://github.com/sympy/sympy/issues
.. [2] https://clang.llvm.org/docs/
.. [3] https://clang.llvm.org/docs/IntroductionToTheClangAST.html

"""

if cin:
    from sympy.codegen.ast import (Variable, IntBaseType, FloatBaseType, String,
        Integer, Float, FunctionPrototype, FunctionDefinition, FunctionCall,
        none, Return, Assignment, Type)
    from sympy.codegen.cnodes import (PreDecrement, PostDecrement,
        PreIncrement, PostIncrement)
    from sympy.core import Add, Mod, Mul, Pow, Rel
    from sympy.logic.boolalg import And, as_Boolean, Not, Or
    from sympy import Symbol, sympify, true, false
    import sys
    import tempfile

    class BaseParser(object):
        """Base Class for the C parser"""

        def __init__(self):
            """Initializes the Base parser creating a Clang AST index"""
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

        The converter object takes the C source code or file as input and
        converts them to SymPy Expressions.
        """

        def __init__(self):
            """Initializes the code converter"""
            super(CCodeConverter, self).__init__()
            self._py_nodes = []

        def parse(self, filenames, flags):
            """Function to parse a file with C source code

            It takes the filename as an attribute and creates a Clang AST
            Translation Unit parsing the file.
            Then the transformation function is called on the transaltion unit,
            whose reults are collected into a list which is returned by the
            function.

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

        def parse_str(self, source, flags):
            """Function to parse a string with C source code

            It takes the source code as an attribute, stores it in a temporary
            file and creates a Clang AST Translation Unit parsing the file.
            Then the transformation function is called on the transaltion unit,
            whose reults are collected into a list which is returned by the
            function.

            Parameters
            ==========

            source : string
                Path to the C file to be parsed

            flags: list
                Arguments to be passed to Clang while parsing the C code

            Returns
            =======

            py_nodes: list
                A list of sympy AST nodes

            """
            file = tempfile.NamedTemporaryFile(mode = 'w+', suffix = '.cpp')
            file.write(source)
            file.seek(0)
            self.tu = self.index.parse(
                file.name,
                args=flags,
                options=cin.TranslationUnit.PARSE_DETAILED_PROCESSING_RECORD
            )
            file.close()
            for child in self.tu.cursor.get_children():
                if child.kind == cin.CursorKind.VAR_DECL:
                    self._py_nodes.append(self.transform(child))
                elif (child.kind == cin.CursorKind.FUNCTION_DECL):
                    self._py_nodes.append(self.transform(child))
                else:
                    pass
            return self._py_nodes

        def transform(self, node):
            """Transformation Function for Clang AST nodes

            It determines the kind of node and calls the respective
            transformation function for that node.

            Raises
            ======

            NotImplementedError : if the transformation for the provided node
            is not implemented

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

            Used to create nodes for variable declarations and assignments with
            values or function call for the respective nodes in the clang AST

            Returns
            =======

            A variable node as Declaration, with the given value or 0 if the
            value is not provided

            Raises
            ======

            NotImplementedError : if called for data types not currently
            implemented

            Notes
            =====

            This function currently only supports basic Integer and Float data
            types

            """
            try:
                children = node.get_children()
                child = next(children)
                #ignoring namespace and type details for the variable
                while child.kind == cin.CursorKind.NAMESPACE_REF:
                    child = next(children)

                while child.kind == cin.CursorKind.TYPE_REF:
                    child = next(children)

                val = self.transform(child)
                # List in case of variable assignment,
                # FunctionCall node in case of a function call
                if (child.kind == cin.CursorKind.INTEGER_LITERAL
                    or child.kind == cin.CursorKind.UNEXPOSED_EXPR):
                    if (node.type.kind == cin.TypeKind.INT):
                        type = IntBaseType(String('integer'))
                        # when only one decl_ref_expr is assigned
                        # e.g., int b = a;
                        if isinstance(val, str):
                            value = Symbol(val)
                        # e.g., int b = true;
                        elif isinstance(val, bool):
                            value = Integer(0) if val == False else Integer(1)
                        # when val is integer or character literal
                        # e.g., int b = 1; or int b = 'a';
                        elif isinstance(val, (Integer, int, Float, float)):
                            value = Integer(val)
                        # when val is combination of both of the above
                        # but in total only two nodes on rhs
                        # e.g., int b = a * 1;
                        else:
                            value = val

                    elif (node.type.kind == cin.TypeKind.FLOAT):
                        type = FloatBaseType(String('real'))
                        # e.g., float b = a;
                        if isinstance(val, str):
                            value = Symbol(val)
                        # e.g., float b = true;
                        elif isinstance(val, bool):
                            value = Float(0.0) if val == False else Float(1.0)
                        # e.g., float b = 1.0;
                        elif isinstance(val, (Integer, int, Float, float)):
                            value = Float(val)
                        # e.g., float b = a * 1.0;
                        else:
                            value = val

                    elif (node.type.kind == cin.TypeKind.BOOL):
                        type = Type(String('bool'))
                        # e.g., bool b = a;
                        if isinstance(val, str):
                            value = Symbol(val)
                        # e.g., bool b = 1;
                        elif isinstance(val, (Integer, int, Float, float)):
                            value = sympify(bool(val))
                        # e.g., bool b = a * 1;
                        else:
                            value = val

                    else:
                        raise NotImplementedError("Only bool, int "
                            "and float are supported")

                elif (child.kind == cin.CursorKind.CALL_EXPR):
                    return Variable(
                        node.spelling
                    ).as_Declaration(
                        value = val
                    )

                # when val is combination of more than two expr and
                # integer(or float) i.e. it has binary operator
                # Or var decl has parenthesis in the rhs(as a parent Clang node)
                # Or it has unary operator
                # Or it has boolean literal on rhs i.e. true or false
                elif (child.kind == cin.CursorKind.BINARY_OPERATOR
                    or child.kind == cin.CursorKind.PAREN_EXPR
                    or child.kind == cin.CursorKind.UNARY_OPERATOR
                    or child.kind == cin.CursorKind.CXX_BOOL_LITERAL_EXPR):
                    if (node.type.kind == cin.TypeKind.INT):
                        type = IntBaseType(String('integer'))
                    elif (node.type.kind == cin.TypeKind.FLOAT):
                        type = FloatBaseType(String('real'))
                    elif (node.type.kind == cin.TypeKind.BOOL):
                        type = Type(String('bool'))
                    else:
                        raise NotImplementedError("Only bool, int "
                            "and float are supported")
                    value = val

                else:
                    raise NotImplementedError("Given "
                        "variable declaration \"{}\" "
                        "is not possible to parse yet!"
                        .format(" ".join(
                            t.spelling for t in node.get_tokens()
                            )
                        ))

            except StopIteration:

                if (node.type.kind == cin.TypeKind.INT):
                    type = IntBaseType(String('integer'))
                    value = Integer(0)
                elif (node.type.kind == cin.TypeKind.FLOAT):
                    type = FloatBaseType(String('real'))
                    value = Float(0.0)
                elif (node.type.kind == cin.TypeKind.BOOL):
                    type = Type(String('bool'))
                    value = false
                else:
                    raise NotImplementedError("Only bool, int "
                            "and float are supported")

            return Variable(
                node.spelling
            ).as_Declaration(
                type = type,
                value = value
            )

        def transform_function_decl(self, node):
            """Transformation Function For Function Declaration

            Used to create nodes for function declarations and definitions for
            the respective nodes in the clang AST

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
                ret_type = IntBaseType(String('integer'))
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

            Used to create parameter nodes for the required functions for the
            respective nodes in the clang AST

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

                # Any namespace nodes can be ignored
                while child.kind in [cin.CursorKind.NAMESPACE_REF,
                                     cin.CursorKind.TYPE_REF,
                                     cin.CursorKind.TEMPLATE_REF]:
                    child = next(children)

                # If there is a child, it is the default value of the parameter.
                lit = self.transform(child)
                if (node.type.kind == cin.TypeKind.INT):
                    val = Integer(lit)
                elif (node.type.kind == cin.TypeKind.FLOAT):
                    val = Float(lit)

                param = Variable(
                    node.spelling
                ).as_Declaration(
                    type = type,
                    value = val
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
            try:
                value = next(node.get_tokens()).spelling
            except StopIteration:
                # No tokens
                value = node.literal
            return int(value)

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
            try:
                value = next(node.get_tokens()).spelling
            except (StopIteration, ValueError):
                # No tokens
                value = node.literal
            return float(value)

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
            """Transformation function for character literal

            Used to get the value of the given character literal.

            Returns
            =======

            val : int
                val contains the ascii value of the character literal

            Notes
            =====

            Only for cases where character is assigned to a integer value,
            since character literal is not in sympy AST

            """
            try:
               value = next(node.get_tokens()).spelling
            except (StopIteration, ValueError):
                # No tokens
               value = node.literal
            return ord(str(value[1]))

        def transform_cxx_bool_literal_expr(self, node):
            """Transformation function for boolean literal

            Used to get the value of the given boolean literal.

            Returns
            =======

            value : bool
                value contains the boolean value of the variable

            """
            try:
                value = next(node.get_tokens()).spelling
            except (StopIteration, ValueError):
                value = node.literal
            return True if value == 'true' else False

        def transform_unexposed_decl(self,node):
            """Transformation function for unexposed declarations"""
            pass

        def transform_unexposed_expr(self, node):
            """Transformation function for unexposed expression

            Unexposed expressions are used to wrap float, double literals and
            expressions

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

            Used to create function call nodes for the function calls present
            in the C code

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
                        param.append(Integer(arg))
                    elif (child.kind == cin.CursorKind.FLOATING_LITERAL):
                        param.append(Float(arg))
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

            These statements are used to wrap different kinds of declararions
            like variable or function declaration
            The function calls the transformer function for the child of the
            given node

            Returns
            =======

            statement : Codegen AST Node
                contains the node returned by the children node for the type of
                declaration

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

        def transform_paren_expr(self, node):
            """Transformation function for Parenthesized expressions

            Returns the result from its children nodes

            """
            return self.transform(next(node.get_children()))

        def transform_unary_operator(self, node):
            """Transformation function for handling unary operators

            Returns
            =======

            unary_expression: Codegen AST node
                    simplified unary expression represented as Codegen AST

            Raises
            ======

            NotImplementedError
                If dereferencing operator(*), address operator(&) or
                bitwise NOT operator(~) is encountered

            """
            # supported operators list
            operators_list = ['+', '-', '++', '--', '!']
            tokens = [token for token in node.get_tokens()]

            # it can be either pre increment/decrement or any other operator from the list
            if tokens[0].spelling in operators_list:
                child = self.transform(next(node.get_children()))
                # (decl_ref) e.g.; int a = ++b; or simply ++b;
                if isinstance(child, str):
                    if tokens[0].spelling == '+':
                        return Symbol(child)
                    if tokens[0].spelling == '-':
                        return Mul(Symbol(child), -1)
                    if tokens[0].spelling == '++':
                        return PreIncrement(Symbol(child))
                    if tokens[0].spelling == '--':
                        return PreDecrement(Symbol(child))
                    if tokens[0].spelling == '!':
                        return Not(Symbol(child))
                # e.g.; int a = -1; or int b = -(1 + 2);
                else:
                    if tokens[0].spelling == '+':
                        return child
                    if tokens[0].spelling == '-':
                        return Mul(child, -1)
                    if tokens[0].spelling == '!':
                        return Not(sympify(bool(child)))

            # it can be either post increment/decrement
            # since variable name is obtained in token[0].spelling
            elif tokens[1].spelling in ['++', '--']:
                child = self.transform(next(node.get_children()))
                if tokens[1].spelling == '++':
                    return PostIncrement(Symbol(child))
                if tokens[1].spelling == '--':
                    return PostDecrement(Symbol(child))
            else:
                raise NotImplementedError("Dereferencing operator, "
                    "Address operator and bitwise NOT operator "
                    "have not been implemented yet!")

        def transform_binary_operator(self, node):
            """Transformation function for handling binary operators

            Returns
            =======

            binary_expression: Codegen AST node
                    simplified binary expression represented as Codegen AST

            Raises
            ======

            NotImplementedError
                If a bitwise operator or
                unary operator(which is a child of any binary
                operator in Clang AST) is encountered

            """
            # get all the tokens of assignment
            # and store it in the tokens list
            tokens = [token for token in node.get_tokens()]

            # supported operators list
            operators_list = ['+', '-', '*', '/', '%','=',
            '>', '>=', '<', '<=', '==', '!=', '&&', '||']

            # this stack will contain variable content
            # and type of variable in the rhs
            combined_variables_stack = []

            # this stack will contain operators
            # to be processed in the rhs
            operators_stack = []

            # iterate through every token
            for token in tokens:
                # token is either '(', ')' or
                # any of the supported operators from the operator list
                if token.kind == cin.TokenKind.PUNCTUATION:

                    # push '(' to the operators stack
                    if token.spelling == '(':
                        operators_stack.append('(')

                    elif token.spelling == ')':
                        # keep adding the expression to the
                        # combined variables stack unless
                        # '(' is found
                        while (operators_stack
                            and operators_stack[-1] != '('):
                            if len(combined_variables_stack) < 2:
                                raise NotImplementedError(
                                    "Unary operators as a part of "
                                    "binary operators is not "
                                    "supported yet!")
                            rhs = combined_variables_stack.pop()
                            lhs = combined_variables_stack.pop()
                            operator = operators_stack.pop()
                            combined_variables_stack.append(
                                self.perform_operation(
                                lhs, rhs, operator))

                        # pop '('
                        operators_stack.pop()

                    # token is an operator (supported)
                    elif token.spelling in operators_list:
                        while (operators_stack
                            and self.priority_of(token.spelling)
                            <= self.priority_of(
                            operators_stack[-1])):
                            if len(combined_variables_stack) < 2:
                                raise NotImplementedError(
                                    "Unary operators as a part of "
                                    "binary operators is not "
                                    "supported yet!")
                            rhs = combined_variables_stack.pop()
                            lhs = combined_variables_stack.pop()
                            operator = operators_stack.pop()
                            combined_variables_stack.append(
                                self.perform_operation(
                                lhs, rhs, operator))

                        # push current operator
                        operators_stack.append(token.spelling)

                    # token is a bitwise operator
                    elif token.spelling in ['&', '|', '^', '<<', '>>']:
                        raise NotImplementedError(
                            "Bitwise operator has not been "
                            "implemented yet!")
                    else:
                        raise NotImplementedError(
                            "Given token {} is not implemented yet!"
                            .format(token.spelling))

                # token is an identifier(variable)
                elif token.kind == cin.TokenKind.IDENTIFIER:
                    combined_variables_stack.append(
                        [token.spelling, 'identifier'])

                # token is a literal
                elif token.kind == cin.TokenKind.LITERAL:
                    combined_variables_stack.append(
                        [token.spelling, 'literal'])

                # token is a keyword, either true or false
                elif (token.kind == cin.TokenKind.KEYWORD
                    and token.spelling in ['true', 'false']):
                    combined_variables_stack.append(
                        [token.spelling, 'boolean'])
                else:
                    raise NotImplementedError(
                        "Given token {} is not implemented yet!"
                        .format(token.spelling))

            # process remaining operators
            while operators_stack:
                if len(combined_variables_stack) < 2:
                    raise NotImplementedError(
                        "Unary operators as a part of "
                        "binary operators is not "
                        "supported yet!")
                rhs = combined_variables_stack.pop()
                lhs = combined_variables_stack.pop()
                operator = operators_stack.pop()
                combined_variables_stack.append(
                    self.perform_operation(lhs, rhs, operator))

            return combined_variables_stack[-1][0]

        def priority_of(self, op):
            """To get the priority of given operator"""
            if op == '=':
                return 1
            if op in ['&&', '||']:
                return 2
            if op in ['<', '<=', '>', '>=', '==', '!=']:
                return 3
            if op in ['+', '-']:
                return 4
            if op in ['*', '/', '%']:
                return 5
            return 0

        def perform_operation(self, lhs, rhs, op):
            """Performs operation supported by the sympy core

            Returns
            =======

            combined_variable: list
                contains variable content and type of variable

            """
            lhs_value = self.get_expr_for_operand(lhs)
            rhs_value = self.get_expr_for_operand(rhs)
            if op == '+':
                return [Add(lhs_value, rhs_value), 'expr']
            if op == '-':
                return [Add(lhs_value, -rhs_value), 'expr']
            if op == '*':
                return [Mul(lhs_value, rhs_value), 'expr']
            if op == '/':
                return [Mul(lhs_value, Pow(rhs_value, Integer(-1))), 'expr']
            if op == '%':
                return [Mod(lhs_value, rhs_value), 'expr']
            if op in ['<', '<=', '>', '>=', '==', '!=']:
                return [Rel(lhs_value, rhs_value, op), 'expr']
            if op == '&&':
                return [And(as_Boolean(lhs_value), as_Boolean(rhs_value)), 'expr']
            if op == '||':
                return [Or(as_Boolean(lhs_value), as_Boolean(rhs_value)), 'expr']
            if op == '=':
                return [Assignment(Variable(lhs_value), rhs_value), 'expr']

        def get_expr_for_operand(self, combined_variable):
            """Gives out SymPy Codegen AST node

            AST node returned is corresponding to
            combined variable passed.Combined variable contains
            variable content and type of variable

            """
            if combined_variable[1] == 'identifier':
                return Symbol(combined_variable[0])
            if combined_variable[1] == 'literal':
                if '.' in combined_variable[0]:
                    return Float(float(combined_variable[0]))
                else:
                    return Integer(int(combined_variable[0]))
            if combined_variable[1] == 'expr':
                return combined_variable[0]
            if combined_variable[1] == 'boolean':
                    return true if combined_variable[0] == 'true' else false

else:
    class CCodeConverter():  # type: ignore
        def __init__(self, *args, **kwargs):
            raise ImportError("Module not Installed")


def parse_c(source):
    """Function for converting a C source code

    The function reads the source code present in the given file and parses it
    to give out SymPy Expressions

    Returns
    =======

    src : list
        List of Python expression strings

    """
    converter = CCodeConverter()
    if os.path.exists(source):
        src = converter.parse(source, flags = [])
    else:
        src = converter.parse_str(source, flags = [])
    return src
