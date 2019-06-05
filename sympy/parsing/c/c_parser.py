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
    __import__kwargs = {
        'fromlist':['cindex']
    }
)

class BaseParser(object):
    def __init__(self):
        self.index = cin.Index.create()

    def diagnostics(self, out):
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
    def __init__(self, name):
        super(CCodeConverter, self).__init__()
        self._py_nodes = []

    def parse(self, filenames, flags):
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

    def parse_text(self, content, flags):
        for f, c in content:
            filename = os.path.abspath(f)

            self.tu = self.index.parse(
                filename,
                args=flags,
                unsaved_files=[(f, c)],
                options=TranslationUnit.PARSE_DETAILED_PROCESSING_RECORD
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
        type = IntBaseType(String('integer'))
        try:
            value = next(node.get_tokens()).spelling
        except StopIteration:
            # No tokens
            value = Integer(node.literal)
        val = [type, value]
        return val

    def transform_floating_literal(self, node):
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
        pass

    def transform_unexposed_expr(self, node):
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
        return node.spelling

    def transform_call_expr(self, node):
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
        return Return(next(node.get_children()).spelling)

    def transform_compound_stmt(self, node):
        try:
            expr = []
            children = node.get_children()
            for child in children:
                expr.append(self.transform(child))
        except StopIteration:
            return None
        return expr

    def transform_decl_stmt(self, node):
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
