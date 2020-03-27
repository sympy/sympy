from sympy.external import import_module
import os

pycp = import_module('pycparser')

"""
This module contains all the necessary Classes and Function used to Parse
C code into SymPy expression
The module serves as a backend for SymPyExpression to parse C code
It is also dependent on pycparser's AST and SymPy's Codegen AST.
The module only supports the features currently supported by the pycparser and
codegen AST which will be updated as the development of codegen AST and this
module progresses.
You might find unexpected bugs and exceptions while using the module, feel free
to report them to the SymPy Issue Tracker

Features Supported
==================

- Variable Declarations (integers and reals)

Notes
=====

The module is dependent on an external dependency which needs to be installed
to use the features of this module.

PyCParser: The C parser which is used to extract an AST from the provided
C source code.

Reference
=========

.. [1] https://github.com/sympy/sympy/issues
.. [2] https://github.com/eliben/pycparser

"""

if pycp:
	from sympy.codegen.ast import (Variable, IntBaseType, FloatBaseType, String,
        Integer, Float)
	
	parser = pycp.CParser()
	class CCodeConverter:

		def __init__(self):
			self.source_string = None
			self.source_file = None
			self.expr_nodes = []

		def parse(self, ast):
			for child in ast.children():
				if isinstance(child[1], pycp.c_ast.Decl):
					self.expr_nodes.append(self.transform_decl(child[1]))

			return self.expr_nodes

		def parse_str(self, c_str):
			self.source_string = c_str
			try:
				ast = parser.parse(c_str)
				return self.parse(ast)
			except pycp.c_parser.ParseError as pe:
				print("Parse error"+str(pe))

		def parse_file(self, c_file):
			self.source_file = c_file
			ast = pycp.c_parser.parse_file(c_file)
			return self.parse(ast)
			

		def transform_decl(self, node):
            # when variable is not initialized
			decl_type = node.type.type.names[0]
			value = None
			if node.init == None:
				if decl_type == 'int':
					type = IntBaseType(String('integer'))
					value = Integer(0)
				elif decl_type == 'float':
					type = FloatBaseType(String('real'))
					value = Float(0.0)
				else:
					raise NotImplementedError("Only integer and float \
						declarations are accepted!")
			else:
				init_type = node.init.type
				print(init_type)
				# when character is assigned to integer variable
				if decl_type == 'int' and init_type == 'char':
					type = IntBaseType(String('integer'))
					value = Integer(self.transform_constant(node))
				if decl_type == 'int':
					type = IntBaseType(String('integer'))
					value = Integer(self.transform_constant(node))
				elif decl_type == 'float':
					type = FloatBaseType(String('real'))
					value = Float(self.transform_constant(node))
				else:
					raise NotImplementedError("Only integer and float \
						declarations are accepted!")
			return Variable(
				node.name
			).as_Declaration(
				type = type,
				value = value
			)

		def transform_constant(self, node):
			if(node.init.type == 'int'):
				return int(node.init.value)
			if(node.init.type == 'float'):
				return float(node.init.value)
			if(node.init.type == 'char'):
				return ord(node.init.value[1])
else:
	class CCodeConverter:
		def __init__(self, *args, **kwargs):
            raise ImportError("Module not Installed")

def parse_c(c_source):
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
    	src = converter.parse_file(source, flags = [])
    else:
    	src = converter.parse_str(source, flags = [])
    return src
