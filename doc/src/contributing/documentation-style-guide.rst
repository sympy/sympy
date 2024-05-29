============================
Documentation Style Guide
============================

General Guidelines
--------------------

Documentation is one of the most highly valued aspects of an open source
project. Documentation teaches users and contributors how to use a project, how
to contribute, and the standards of conduct within an open source community.
But according to GitHub’s `Open Source Survey
<https://opensourcesurvey.org/2017/>`_, incomplete or confusing documentation is
the most commonly encountered problem in open source. This style guide aims to
change that.

The purpose of this style guide is to provide the SymPy community with a set of
style and formatting guidelines that can be utilized and followed when writing
SymPy documentation. Adhering to the guidelines offered in this style guide
will bring greater consistency and clarity to SymPy’s documentation, supporting
its mission to become a full-featured, open source computer algebra system
(CAS).

The SymPy documentation found at `docs.sympy.org
<https://docs.sympy.org/latest/index.html>`_ is generated from docstrings in the
source code and dedicated narrative documentation files in the `doc/src
directory <https://github.com/sympy/sympy/tree/master/doc/src>`_. Both are
written in `reStructuredText <https://docutils.sourceforge.io/rst.html>`_ format
extended by `Sphinx <https://www.sphinx-doc.org/en/master/>`_.

The documentation contained in the `doc/src directory
<https://github.com/sympy/sympy/tree/master/doc/src>`_ and the docstrings
embedded in the Python source code are processed by Sphinx and various Sphinx
extensions. This means that the documentation source format is specified by the
documentation processing tools. The SymPy Documentation Style Guide provides
both the essential elements for writing SymPy documentation as well as any
deviations in style we specify relative to these documentation processing tools.
The following lists the processing tools:

* reStructuredText: Narrative documentation files and documentation strings
  embedded in Python code follow the reStructuredText format. Advanced features
  not described in this document can be found at
  https://docutils.sourceforge.io/rst.html.

* Sphinx: Sphinx includes additional default features for the
  reStructuredText specification that are described at: https://www.sphinx-doc.org/en/master.

* Sphinx extensions included with Sphinx that we enable:

  * ``sphinx.ext.autodoc``: Processes Python source code files for the
    associated documentation strings to automatically generate pages containing
    the Application Programming Interface (API). See section on calling autodoc
    directives in this document to get started. More information is at:
    https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html.
  * ``sphinx.ext.graphviz``: Provides a directive for adding Graphviz graphics.
    See https://www.sphinx-doc.org/en/master/usage/extensions/graphviz.html for
    more info.
  * ``sphinx.ext.mathjax``: Causes math written in LaTeX to display using
    MathJax in the HTML version of the documentation. More information is at:
    https://www.sphinx-doc.org/en/master/usage/extensions/math.html#module-sphinx.ext.mathjax.
    *No bearing on documentation source format.*
  * ``sphinx.ext.linkcode``: Causes links to source code to direct to the
    related files on Github. More information is at:
    https://www.sphinx-doc.org/en/master/usage/extensions/linkcode.html. *No
    bearing on documentation source format.*

* Sphinx extensions that are not included with Sphinx that we enable:

  * ``numpydoc``: Processes docstrings written in the "numpydoc" format, see
    https://numpydoc.readthedocs.io/en/stable/. We recommend the subset of numpydoc
    formatting features in this document. (Note that we currently use an older
    modified fork of numpydoc, which is included in the SymPy source code.)
  * ``sphinx_math_dollar``: Allows math to be delimited with dollar signs
    instead of reStructuredText directives (e.g., ``$a^2$`` instead of
    ``:math:`a^2```). See https://www.sympy.org/sphinx-math-dollar/ for more info.
  * ``matplotlib.sphinxext.plot_directive``: Provides directives for included
    matplotlib generated figures in reStructuredText. See
    https://matplotlib.org/devel/plot_directive.html for more info.

Everything supported by the above processing tools is available for use in the
SymPy documentation, but this style guide supersedes any recommendations made
in the above documents. Note that we do not follow PEP 257 or the
www.python.org documentation recommendations.

If you are contributing to SymPy for the first time, please read our
:doc:`introduction-to-contributing` page as
well as this guide.

Types of Documentation
------------------------

There are four main locations where SymPy’s documentation can be found:

**SymPy Website** https://www.sympy.org/

The SymPy website’s primary function is to advertise the software to users and
developers. It also serves as an initial location to point viewers to other
relevant resources on the web. The SymPy website has basic information on SymPy
and how to obtain it, as well as examples to advertise it to users, but it does
not have technical documentation. The source files are located in the SymPy
`webpage directory <https://github.com/sympy/sympy.github.com>`_. Appropriate
items for the website are:

* General descriptions of what SymPy and the SymPy community are
* Explanations/demonstrations of major software features
* Listings of other major software that uses SymPy
* Getting started info for users (download and install instructions)
* Getting started info for developers
* Where users can get help and support on using SymPy
* News about SymPy

**SymPy Documentation** https://docs.sympy.org

This is the main place where users go to learn how to use SymPy. It contains a
tutorial for SymPy as well as technical documentation for all of the modules.
The source files are hosted in the main SymPy repository in the `doc directory
<https://github.com/sympy/sympy/tree/master/doc>`_ at and are built using the
`Sphinx site generator <https://www.sphinx-doc.org/en/master/>`_ and uploaded
to the docs.sympy.org site automatically. There are two primary types of pages
that are generated from different source files in the docs directory:

* Narrative Pages: reStructuredText files that correspond to manually written
  documentation pages not present in the Python source code. Examples are the
  `tutorial RST files
  <https://github.com/sympy/sympy/tree/master/doc/src/tutorials>`_. In general,
  if your documentation is not API documentation it belongs in a narrative page.
* API Documentation Pages: reStructuredText files that contain directives that
  generate the Application Programming Interface documentation. These are
  automatically generated from the SymPy Python source code.

**SymPy Source Code** https://github.com/sympy/sympy

Most functions and classes have documentation written inside it in the form of a
docstring, which explains the function and includes examples called doctests.
The purpose of these docstrings are to explain the API of that class or
function. The doctests examples are tested as part of the test suite, so that we
know that they always produce the output that they say that they do. Here is an
`example docstring
<https://github.com/sympy/sympy/blob/b176f6a1d9890b42dc361857c887992315e3d5ad/sympy/functions/elementary/complexes.py#L22-L47>`_.
Most docstrings are also automatically included in the Sphinx documentation
above, so that they appear on the SymPy Documentation website. Here is that
:obj:`same docstring <.im>` on the SymPy website. The docstrings are formatted
in a specific way so that Sphinx can render them correctly for the docs
website. The SymPy sources all contain sparse technical documentation in the
form of source code comments, although this does not generally constitute
anything substantial and is not displayed on the documentation website.

**SymPy Wiki** https://github.com/sympy/sympy/wiki

The SymPy Wiki can be edited by anyone without review. It contains various
types of documentation, including:

* High-level developer documentation (for example: https://github.com/sympy/sympy/wiki/Args-Invariant)
* Release notes (for example: https://github.com/sympy/sympy/wiki/Release-Notes-for-1.5)
* Various pages that different contributors have added

Narrative Documentation Guidelines
-----------------------------------

Extensive documentation, or documentation that is not centered around an API
reference, should be written as a narrative document in the Sphinx docs (located
in the `doc/src directory
<https://github.com/sympy/sympy/tree/master/doc/src>`_). The narrative documents
do not reside in the Python source files, but as standalone restructured files
in the doc/src directory. SymPy’s narrative documentation is defined as the
collective documents, tutorials, and guides that teach users how to use SymPy.
Reference documentation should go in the docstrings and be pulled into the RST
with autodoc. The RST itself should only have narrative style documentation
that is not a reference for a single specific function.

Documentation using Markdown
----------------------------

Narrative documentation can be written using either Restructured Text
(``.rst``) or Markdown (``.md``). Markdown documentation uses `MyST
<https://myst-parser.readthedocs.io/en/latest/index.html>`_. See `this guide
<https://myst-parser.readthedocs.io/en/latest/syntax/syntax.html>`_ for more
information on how to write documents in Markdown. Markdown is only supported
for narrative documentation. Docstrings should continue to use RST syntax. Any
part of this style guide that is not specific to RST syntax should still apply
to Markdown documents.


.. _style_guide_best_practices_for_writing_documentation:

Best Practices for Writing Documentation
----------------------------------------

Please follow these formatting, style, and tone preferences when writing
documentation.

Formatting Preferences
^^^^^^^^^^^^^^^^^^^^^^

In order for math and code to render correctly on the SymPy website, please
follow these formatting guidelines.

.. _style_guide_math_formatting:

Math
~~~~

Text that is surrounded by dollar signs $ _ $ will be rendered as LaTeX math.
Any text that is meant to appear as LaTeX math should be written as ``$math$``.
In the HTML version of the docs, MathJax will render the math.

**Example**

::

    The Bessel $J$ function of order $\nu$ is defined to be the function
    satisfying Bessel’s differential equation.

.. _style_guide_latex_recommendations:

LaTeX Recommendations
~~~~~~~~~~~~~~~~~~~~~

* If a docstring has any LaTeX, be sure to make it "raw." See the
  :ref:`Docstring Formatting <style_guide_docstring_formatting>` section for
  details.
* If you are not sure how to render something, you can use the SymPy
  :func:`~.latex` function. But be sure to strip out the unimportant parts (the
  bullet points below).
* Avoid unnecessary ``\left`` and ``\right`` (but be sure to use them when they
  are required).
* Avoid unnecessary ``{}``. (For example, write ``x^2`` instead of ``x^{2}``.)
* Use whitespace in a way that makes the equation easiest to read.
* Always check the final rendering to make sure it looks the way you expect it
  to.
* The HTML documentation build will not fail if there is invalid math, but
  rather it will show as an error on the page. However, the PDF build, which
  is run on GitHub Actions on pull requests, will fail. If the LaTeX PDF build
  fails on CI, there is likely an issue with LaTeX math somewhere.

**Examples**

Correct::

    \int \sin(x)\,dx

Incorrect::

    \int \sin{\left( x\right)}\, dx

For more in-depth resources on how to write math in LaTeX, see:

* https://math.meta.stackexchange.com/questions/5020/mathjax-basic-tutorial-and-quick-reference
* https://en.wikibooks.org/wiki/LaTeX/Mathematics
* https://www.overleaf.com/learn/latex/Mathematical_expressions

Code
~~~~

Text that should be printed verbatim, such as code, should be surrounded by a
set of double backticks ``like this``.

**Example**

::

    To use this class, define the ``_rewrite()`` and ``_expand()`` methods.

Sometimes a variable will be the same in both math and code, and can even
appear in the same paragraph, making it difficult to know if it should be
formatted as math or code. If the sentence in question is discussing
mathematics, then LaTeX should be used, but if the sentence is discussing the
SymPy implementation specifically, then code should be used.

In general, the rule of thumb is to consider if the variable in question were
something that rendered differently in code and in math. For example, the Greek
letter α would be written as ``alpha`` in code and ``$\alpha$`` in LaTeX. The
reason being that ``$\alpha$`` cannot be used in contexts referring to Python
code because it is not valid Python, and conversely ``alpha`` would be
incorrect in a math context because it does not render as the Greek letter (α).

**Example**

::

    class loggamma(Function):
        r"""
        The ``loggamma`` function implements the logarithm of the gamma
        function (i.e, $\log\Gamma(x)$).

        """

Variables listed in the parameters after the function name should, in written
text, be italicized using Sphinx emphasis with asterisks like ``*this*``.

**Example**

::

    def stirling(n, k, d=None, kind=2, signed=False):
        """
        ...

        The first kind of Stirling number counts the number of permutations of
        *n* distinct items that have *k* cycles; the second kind counts the
        ways in which *n* distinct items can be partitioned into *k* parts.
        If *d* is given, the "reduced Stirling number of the second kind" is
        returned: $S^{d}(n, k) = S(n - d + 1, k - d + 1)$ with $n \ge k \ge d$.
        This counts the ways to partition $n$ consecutive integers into $k$
        groups with no pairwise difference less than $d$.

        """

Note that in the above example, the first instances of *n* and *k* are
referring to the input parameters of the function ``stirling``. Because they
are Python variables but also parameters listed by themselves, they are
formatted as parameters in italics. The last instances of $n$ and $k$ are
talking about mathematical expressions, so they are formatted as math.

If a variable is code, but is also a parameter written by itself, the parameter
formatting takes precedence and it should be rendered in italics. However, if a
parameter appears in a larger code expression it should be within double
backticks to be rendered as code. If a variable is only code and not a
parameter as well, it should be in double backticks to be rendered as code.

Please note that references to other functions in SymPy are handled differently
from parameters or code. If something is referencing another function in SymPy,
the cross-reference reStructuredText syntax should be used. See the section on
:ref:`Cross-Referencing <style_guide_cross-referencing>` for more information.

Headings
~~~~~~~~

Section headings in reStructuredText files are created by underlining (and
optionally overlining) the section title with a punctuation character at least
as long as the text.

Normally, there are no heading levels assigned to certain characters as the
structure is determined from the succession of headings. However, for SymPy's
documentation, here is a suggested convention:

``===`` with overline: title (top level heading)

``===`` heading 1

``---`` heading 2

``^^^`` heading 3

``~~~`` heading 4

``"""`` heading 5

Style Preferences
^^^^^^^^^^^^^^^^^

Spelling and Punctuation
~~~~~~~~~~~~~~~~~~~~~~~~

All narrative writing in SymPy follows American spelling and punctuation
standards. For example, “color” is preferred over “colour” and commas should be
placed inside of quotation marks.

**Examples**

::

    If the ``line_color`` aesthetic is a function of arity 1, then the coloring
    is a function of the x value of a point.

    The term "unrestricted necklace," or "bracelet," is used to imply an object
    that can be turned over or a sequence that can be reversed.

If there is any ambiguity about the spelling of a word, such as in the case of
a function named after a person, refer to the spelling of the actual SymPy
function.

For example, Chebyshev polynomials are named after Pafnuty Lvovich Tchebychev,
whose name is sometimes transliterated from Russian to be spelled with a “T,”
but in SymPy it should always be spelled “Chebyshev” to refer to the SymPy
function.

**Example**

::

    class chebyshevt(OrthogonalPolynomial):
        r"""
        Chebyshev polynomial of the first kind, $T_n(x)$
        ...

        """

Capitalization
~~~~~~~~~~~~~~

Title case capitalization is preferred in all SymPy headings.

**Example**

::

    What is Symbolic Computation?
    -----------------------------

Tone Preferences
^^^^^^^^^^^^^^^^

Across SymPy documentation, please write in:

* The present tense (e.g., In the following section, we are going to learn...)
* The first-person inclusive plural (e.g., We did this the long way, but now we
  can try it the short way...)
* Use the generic pronoun “you” instead of “one.” Or use “the reader” or “the
  user.” (e.g., You can access this function by... The user can then access
  this function by...)
* Use the gender-neutral pronoun “they” instead of “he” or “she.” (e.g., A good
  docstring tells the user exactly what they need to know.)

Avoid extraneous or belittling words such as “obviously,” “easily,” “simply,”
“just,” or “straightforward.”

Avoid unwelcoming or judgement-based phrases like “That is wrong.” Instead use
friendly and inclusive language like “A common mistake is...”

Avoid extraneous phrases like, “we just have to do one more thing.”
