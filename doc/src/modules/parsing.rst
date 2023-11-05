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

The current implementations are experimental. The behavior, parser backend(s) and
API might change in the future. Unlike some of the other parsers, `\mathrm{\LaTeX}` is
designed as a *type-setting* language, not a *computer algebra system* and so
can contain typographical conventions that might be interpreted multiple ways.

`\mathrm{\LaTeX}` Parsing Functions Reference
---------------------------------------------

.. autofunction:: sympy.parsing.latex.parse_latex

ANTLR Backend
^^^^^^^^^^^^^

The ANTLR-based `\mathrm{\LaTeX}` parser was ported from
`latex2sympy <https://github.com/augustt198/latex2sympy>`_. While functional
and its API should remain stable, the parsing behavior or backend may change in
future releases.

.. _ANTLR parser caveats:

ANTLR `\mathrm{\LaTeX}` Parser Caveats
""""""""""""""""""""""""""""""""""""""

In its current definition, the parser may fail to fully parse
an expression, yet not throw a warning::

    parse_latex(r'x -')

will simply find ``x``. What is covered by this behavior will almost certainly
change between releases, and become stricter, more relaxed, or some mix.

Lark Backend
^^^^^^^^^^^^

The Lark-based LaTeX parser is newer, and is intended to eventually completely
replace the ANTLR-based parser. It has most of the features that the ANTLR-based
parser provides, with some extras.

Lark `\mathrm{\LaTeX}` Parser Features
""""""""""""""""""""""""""""""""""""""

One thing to note is that the Lark backend does not support ill-formed expressions,
and it does not try to fix any sort of common mistakes that may have occured. For
example, as mentioned in :ref:`the earlier section <ANTLR parser caveats>`, the
ANTLR-based parser would simply find ``x`` if we run::

    parse_latex(r'x -', backend='ANTLR')

However, running::

    parse_latex(r'x -', backend='Lark')

will raise an ``lark.exceptions.UnexpectedEOF`` exception.

Apart from that, there are a couple of extra things that the Lark-based parser
supports that the ANTLR-based parser does not. They are:

1. Detecting ambiguous expressions, and
2. Allowing user-customization of the `\mathrm{\LaTeX}` grammar at runtime.

Expressions like `f(x)` are technically ambiguous `\mathrm{\LaTeX}` expressions
because the `f` might be a function or a variable name. Lark has the capability to
point out these ambiguities and notify the user, or even return all possible
interpretations.

The Lark-based parser exposes a number of internals which allow the user to customize
the parser's behavior. For example, the user can specify their own `\mathrm{\LaTeX}`
grammar by passing the path to the grammar file to the ``LarkLaTeXParser`` while
instantiating the parser.

The user can also specify their own custom transformer class to the `LarkLaTeXParser`
class.

The two examples mentioned above can be found in the
`test_custom_latex.py <https://github.com/sympy/sympy/blob/395e820b114d2b169483354f1f4ee2f439faa292/sympy/parsing/tests/test_custom_latex.py>`_
file.

Lark `\mathrm{\LaTeX}` Parser Capabilities
""""""""""""""""""""""""""""""""""""""""""

In order to use the Lark-based LaTeX parser, it is important to know what it can
and cannot do. As the parser is still experimental, it supports many things, but
some features are still only partially implemented, or not available.

As such, we will list the types of expressions that it can parse, and then list
some expression types of interest where it may fail.

Here is a list of the things which are supported:

* Symbols which consist of one letter, e.g., ``a``, ``b``, ``x``, etc.
  Greek symbols and symbols with subscripts are also supported. Numbers are also
  supported, as is ``\infty``.
* Symbols with multiple letters are supported, as long as they are wrapped in
  ``\mathit``.
* Expressions with `+`, `-`, `*`, `/`, and alternative operators like ``\cdot``,
  ``\times``, ``\div``, etc. If two expressions are next to each other, like `xy`
  or `(\sin x)(\cos t)`, then it is treated as implicit multiplication.
* Relations with `<`, `>`, `\le`, `\ge`, `=`, and `\ne`.
* Commonly used functions like

    * Square roots,
    * factorials,
    * complex conjugation (like `\overline{z}`)
    * `\log`,
    * `\ln`,
    * `\exp`,
    * absolute value (e.g., `|x|`). Note that `||x||` is parsed as ``Abs(Abs(x))``.
    * floor (e.g., `\lfloor x \rfloor`) and ceiling (e.g., `\lceil x \rceil`),
    * `\min` and `\max`.

* All the trigonometric functions and their inverses trigonometric functions.
  Powers like ``\sin^4`` are also supported. The power `-1` is interpreted as the
  inverse function (i.e., ``\sin^{-1} x`` is interpreted as ``\arcsin x``).
* Hyperbolic trigonometric functions (currently only `\sinh`, `\cosh`, and
  `\tanh`) and their inverses. As mentioned in the previous point, powers like
  ``\tanh^2`` are also supported, and `-1` is interpreted as the inverse function
  (i.e., ``\tanh^{-1} x`` is interpreted as ``\arctanh x``).
* ``AppliedFunctions``, like `f(x, y, z)`.
* All types of fractions (``\frac``, ``\tfrac``, ``\dfrac``, ``\nicefrac``) and
  binomials (``\binom``, ``\tbinom``, ``\dbinom``) are supported.
* Integrals, both definite and indefinite. When the integrand is a fraction,
  having the differential in the numerator is allowed. The differential is
  allowed to be ``d``, ``\text{d}``, or ``\mathrm{d}``.
* Derivatives in one variable. I.e., things like `\dfrac{d}{dx} (\sin x)`.
  Higher order derivatives and partial derivatives are not supported yet.
* Limits in one variable. E.g., `\lim\limits_{t\to 3^{+}} \sin t`.
* Sums and products with simple conditions. For example, `\sum\limits_{k=0}^n k^2`
  is allowed because the condition on `k` is simple. An expression like
  `\sum\limits_{d \mid n} d^2` is not allowed because the condition on `d` in
  the subscript is complicated. Expressions with the index variable specified
  in the superscript are also allowed. For example, `\prod\limits_{k=0}^{k=n} k^2`
  is parsed correctly.
* Bra (e.g., `| x \rangle`), and Ket (e.g., `\langle p |`) notation. Parsing
  Inner (e.g., `\langle x | y \rangle`) and Outer Products (e.g.,
  `| y \rangle \langle x |`) is also supported.

Here is a(n incomplete) list of things which are currently not supported, which
may be added in the future:

* Matrices. Stuff like ``\begin{env}...\end{env}``, where ``env`` is any of
  ``matrix``, ``bmatrix``, ``pmatrix``, ``smallmatrix``, and ``array``.
* Matrix operations like matrix-matrix addition, scalar-matrix multiplication,
  matrix-matrix multiplication.
* Higher order derivatives and partial derivatives.
* Double and triple integrals.


Lark `\mathrm{\LaTeX}` Parser Functions
"""""""""""""""""""""""""""""""""""""""

.. autofunction:: sympy.parsing.latex.parse_latex_lark

Lark `\mathrm{\LaTeX}` Parser Classes
"""""""""""""""""""""""""""""""""""""

.. autoclass:: sympy.parsing.latex.lark.LarkLaTeXParser

.. autoclass:: sympy.parsing.latex.lark.TransformToSymPyExpr

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
`ANTLR4 <https://www.antlr.org/>`_,
but to use the parser, you only need the ``antlr4`` Python package available.

Depending on your package manager, you can install the right package with, for
example, ``pip``::

    $ pip install antlr4-python3-runtime==4.11

or ``conda``::

    $ conda install -c conda-forge antlr-python-runtime==4.11

The C parser depends on ``clang`` and the Fortran parser depends on ``LFortran``.
You can install these packages using::

    $ conda install -c conda-forge lfortran clang
