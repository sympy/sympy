Parsing input
=============

Parsing Functions Reference
---------------------------

.. autofunction:: sympy.parsing.mathematica.mathematica

.. autofunction:: sympy.parsing.maxima.parse_maxima

.. autofunction:: sympy.parsing.sympy_parser.eval_expr

.. autofunction:: sympy.parsing.sympy_parser.parse_expr

.. autofunction:: sympy.parsing.sympy_parser.stringify_expr

.. autofunction:: sympy.parsing.sympy_tokenize.any

.. autofunction:: sympy.parsing.sympy_tokenize.generate_tokens

.. autofunction:: sympy.parsing.sympy_tokenize.group

.. autofunction:: sympy.parsing.sympy_tokenize.maybe

.. autofunction:: sympy.parsing.sympy_tokenize.printtoken

.. autofunction:: sympy.parsing.sympy_tokenize.tokenize

.. autofunction:: sympy.parsing.sympy_tokenize.untokenize


Parsing Exceptions Reference
----------------------------

.. autoclass:: sympy.parsing.sympy_tokenize.TokenError

.. autoclass:: sympy.parsing.sympy_tokenize.StopTokenizing

Parsing Transformations Reference
---------------------------------

A transformation is a function that accepts the arguments ``tokens,
local_dict, global_dict`` and returns a list of transformed tokens. They can
be used by passing a list of functions to :py:func:`parse_expr` and are
applied in the order given.

.. autofunction:: sympy.parsing.sympy_parser.convert_xor

.. autofunction:: sympy.parsing.sympy_parser.function_exponentiation

.. autofunction:: sympy.parsing.sympy_parser.implicit_application

.. autofunction:: sympy.parsing.sympy_parser.implicit_multiplication

.. autofunction:: sympy.parsing.sympy_parser.implicit_multiplication_application

.. autofunction:: sympy.parsing.sympy_parser.split_symbols

.. autofunction:: sympy.parsing.sympy_parser.split_symbols_custom

.. autodata:: sympy.parsing.sympy_parser.standard_transformations

.. autofunction:: sympy.parsing.sympy_parser.rationalize


These are included in
:data:``sympy.parsing.sympy_parser.standard_transformations`` and generally
don't need to be manually added by the user.

.. autofunction:: sympy.parsing.sympy_parser.auto_number

.. autofunction:: sympy.parsing.sympy_parser.auto_symbol

.. autofunction:: sympy.parsing.sympy_parser.factorial_notation
