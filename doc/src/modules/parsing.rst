=======
Parsing
=======

.. module:: sympy.parsing

Parsing Functions Reference
---------------------------

.. autofunction:: sympy.parsing.sympy_parser.parse_expr

.. autofunction:: sympy.parsing.sympy_parser.stringify_expr

.. autofunction:: sympy.parsing.sympy_parser.eval_expr

.. autofunction:: sympy.parsing.maxima.parse_maxima

.. autofunction:: sympy.parsing.mathematica.parse_mathematica

Parsing Transformations Reference
---------------------------------

A transformation is a function that accepts the arguments ``tokens,
local_dict, global_dict`` and returns a list of transformed tokens. They can
be used by passing a list of functions to :py:func:`~.parse_expr` and are
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
:data:`sympy.parsing.sympy_parser.standard_transformations` and generally
don't need to be manually added by the user.

.. autofunction:: sympy.parsing.sympy_parser.lambda_notation

.. autofunction:: sympy.parsing.sympy_parser.auto_symbol

.. autofunction:: sympy.parsing.sympy_parser.repeated_decimals

.. autofunction:: sympy.parsing.sympy_parser.auto_number

.. autofunction:: sympy.parsing.sympy_parser.factorial_notation

Experimental `\mathrm{\LaTeX}` Parsing
--------------------------------------

`\mathrm{\LaTeX}` parsing was ported from
`latex2sympy <https://github.com/augustt198/latex2sympy>`_. While functional
and its API should remain stable, the parsing behavior or backend may change in
future releases.

`\mathrm{\LaTeX}` Parsing Caveats
---------------------------------

The current implementation is experimental. The behavior, parser backend and
API might change in the future. Unlike some of the other parsers, `\mathrm{\LaTeX}` is
designed as a *type-setting* language, not a *computer algebra system* and so
can contain typographical conventions that might be interpreted multiple ways.

In its current definition, the parser will at times will fail to fully parse
the expression, but not throw a warning::

    parse_latex(r'x -')

Will simply find ``x``. What is covered by this behavior will almost certainly
change between releases, and become stricter, more relaxed, or some mix.


`\mathrm{\LaTeX}` Parsing Functions Reference
---------------------------------------------

.. autofunction:: sympy.parsing.latex.parse_latex

`\mathrm{\LaTeX}` Parsing Exceptions Reference
----------------------------------------------

.. autoclass:: sympy.parsing.latex.LaTeXParsingError
   :members:

SymPy Expression Reference
--------------------------

.. module:: sympy.parsing.sym_expr

.. autoclass:: SymPyExpression
   :members:

Runtime Installation
--------------------

The currently-packaged LaTeX parser backend is partially generated with
`ANTLR4 <http://antlr4.org>`_,
but to use the parser, you only need the ``antlr4`` Python package available.

Depending on your package manager, you can install the right package with, for
example, ``pip``::

    $ pip install antlr4-python3-runtime==4.11

or ``conda``::

    $ conda install -c conda-forge antlr-python-runtime==4.11

The C parser depends on ``clang`` and the Fortran parser depends on ``LFortran``.
You can install these packages using::

    $ conda install -c conda-forge lfortran clang
