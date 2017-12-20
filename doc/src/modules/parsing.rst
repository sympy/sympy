Parsing input
=============

Parsing Functions Reference
---------------------------

.. autofunction:: sympy.parsing.sympy_parser.parse_expr

.. autofunction:: sympy.parsing.sympy_parser.stringify_expr

.. autofunction:: sympy.parsing.sympy_parser.eval_expr

.. autofunction:: sympy.parsing.sympy_tokenize.printtoken

.. autofunction:: sympy.parsing.sympy_tokenize.tokenize

.. autofunction:: sympy.parsing.sympy_tokenize.untokenize

.. autofunction:: sympy.parsing.sympy_tokenize.generate_tokens

.. autofunction:: sympy.parsing.sympy_tokenize.group

.. autofunction:: sympy.parsing.sympy_tokenize.any

.. autofunction:: sympy.parsing.sympy_tokenize.maybe

.. autofunction:: sympy.parsing.maxima.parse_maxima

.. autofunction:: sympy.parsing.mathematica.mathematica

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

.. autodata:: sympy.parsing.sympy_parser.standard_transformations

.. autofunction:: sympy.parsing.sympy_parser.split_symbols

.. autofunction:: sympy.parsing.sympy_parser.split_symbols_custom

.. autofunction:: sympy.parsing.sympy_parser.implicit_multiplication

.. autofunction:: sympy.parsing.sympy_parser.implicit_application

.. autofunction:: sympy.parsing.sympy_parser.function_exponentiation

.. autofunction:: sympy.parsing.sympy_parser.implicit_multiplication_application

.. autofunction:: sympy.parsing.sympy_parser.rationalize

.. autofunction:: sympy.parsing.sympy_parser.convert_xor

These are included in
:data:``sympy.parsing.sympy_parser.standard_transformations`` and generally
don't need to be manually added by the user.

.. autofunction:: sympy.parsing.sympy_parser.factorial_notation

.. autofunction:: sympy.parsing.sympy_parser.auto_symbol

.. autofunction:: sympy.parsing.sympy_parser.auto_number


Experimental `\LaTeX` Parsing
-----------------------------

`\LaTeX` parsing was ported from
`latex2sympy <https://github.com/augustt198/latex2sympy>`_. While functional
and its API should remain stable, the parsing behavior or backend may change in
future releases.

`\LaTeX` Parsing Functions Reference
------------------------------------

.. autofunction:: sympy.parsing.latex.parse_latex

`\LaTeX` Parsing Exceptions Reference
-------------------------------------

.. autoclass:: sympy.parsing.latex.LaTeXSyntaxError

Runtime Installation
--------------------

The currently-packaged backend is partially generated with
`ANTLR4 <http://antlr4.org>`_,
but to use the parser, you only need the ``antlr4`` python package, provided by
one of::

    $ conda install -c conda-forge antlr-python-runtime
    $ pip install antlr4-python3-runtime
    $ pip install antlr4-python2-runtime
