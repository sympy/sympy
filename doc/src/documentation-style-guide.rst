===============================
SymPy Documentation Style Guide
===============================

**A Writing Resource for Documentation and Docstrings**

General Guidelines
==================

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
written in `reStructuredText <http://docutils.sourceforge.net/rst.html>`_ format
extended by `Sphinx <http://www.sphinx-doc.org/en/master/>`_.

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
  http://docutils.sourceforge.net/rst.html.

* Sphinx: Sphinx includes additional default features for the
  reStructuredText specification that are described at: http://www.sphinx-doc.org/.

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
    https://numpydoc.readthedocs.io. We recommend the subset of numpydoc
    formatting features in this document. (Note that we currently use an older
    modified fork of numpydoc, which is included in the SymPy source code.)
  * ``sphinx_math_dollar``: Allows math to be delimited with dollar signs
    instead of reStructuredText directives (e.g., ``$a^$`` instead of
    ``:math:a^2``). See https://www.sympy.org/sphinx-math-dollar/ for more info.
  * ``matplotlib.sphinxext.plot_directive``: Provides directives for included
    matplotlib generated figures in reStructuredText. See
    https://matplotlib.org/devel/plot_directive.html for more info.
  * ``sympylive``: Adds a button on each example in the HTML documentation that
    opens the example in SymPy Live. *No bearing on documentation source
    format.*

Everything supported by the above processing tools is available for use in the
SymPy documentation, but this style guide supersedes any recommendations made
in the above documents. Note that we do not follow PEP 257 or the
www.python.org documentation recommendations.

If you are contributing to SymPy for the first time, please read our
`Introduction to Contributing
<https://github.com/sympy/sympy/wiki/Introduction-to-contributing>`_ page as
well as this guide.

Types of Documentation
======================

There are four main locations where SymPy’s documentation can be found:

**SymPy Website** https://sympy.org

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
`Sphinx site generator <http://www.sphinx-doc.org/en/master/>`_ and uploaded to
the docs.sympy.org site automatically. The docs website also contains a built-
in shell (SymPy Live) that allows users to interactively execute examples.
There are two primary types of pages that are generated from different source
files in the docs directory:

* Narrative Pages: reStructuredText files that correspond to manually written
  documentation pages not present in the Python source code. Examples are the
  `tutorial RST files
  <https://github.com/sympy/sympy/tree/master/doc/src/tutorial>`_. In general,
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
* Guides for new contributors (for example: https://github.com/sympy/sympy/wiki/Introduction-to-contributing)
* Development policies (for example: https://github.com/sympy/sympy/wiki/Python-version-support-policy)
* Release notes (for example: https://github.com/sympy/sympy/wiki/Release-Notes-for-1.5)
* Various pages that different contributors have added

Getting Started
===============

The first step to contributing to the code base is creating your development
environment. Please find instructions on how to create your development
environment in our `Development Workflow – Create Your Environment
<https://github.com/sympy/sympy/wiki/Development-workflow#create-your-environment>`_
guide.

Once you have created your development environment, follow these steps:

1. Installation
---------------

**Debian/Ubuntu**

For Debian/Ubuntu, install the prerequisites::

   apt-get install python-sphinx texlive-latex-recommended dvipng librsvg2-bin
   imagemagick docbook2x graphviz
   python -m pip install sphinx-math-dollar

And do::

   make html

If you get mpmath error, install python-mpmath package::

   apt-get install python-mpmath

If you get matplotlib error, install python-matplotlib package::

   apt-get install python-matplotlib

And to view it, do::

   firefox _build/html/index.html

**Fedora**

For Fedora (and maybe other RPM-based distributions), install the
prerequisites::

   dnf install python3-sphinx librsvg2 ImageMagick docbook2X texlive-dvipng-bin
   texlive-scheme-medium librsvg2-tools python -m pip install sphinx-math-dollar

After that, run::

   make html

If you get mpmath error, install python3-mpmath package::

   dnf install python3-mpmath

If you get matplotlib error, install python3-matplotlib package::

   dnf install python3-matplotlib

And view it at::

   _build/html/index.html

**Mac**

For Mac, first install homebrew: https://brew.sh/

Then install these packages with homebrew::

   brew install imagemagick graphviz docbook librsvg

Install these packages with either pip or conda::

   python -m pip install mpmath matplotlib sphinx sphinx-math-dollar

Or::

   conda install -c conda-forge mpmpath matplotlib sphinx sphinx-math-dollar

**Windows 10**

Making your Sphinx build successful on the Windows system is tricky because
some dependencies like ``dvipng`` or ``docbook2x`` are not available.

For Windows 10, however, the Windows Subsystem for Linux can be a possible
workaround solution, and you can install Ubuntu shell on your Windows system
after following the tutorial below:

https://github.com/MicrosoftDocs/WSL/blob/live/WSL/install-win10.md

In your command prompt, run ``ubuntu`` to transfer to Linux terminal, and
follow the Debian/Ubuntu tutorial above to install the dependencies, and then
you can run ``make html`` to build. (Note that you also have to install
``make`` via ``apt-get install make``.)

If you want to change the directory in your prompt to your working folder of
SymPy in the Windows file system, you can prepend ``cd /mnt/`` to your file
path in Windows, and run in your shell to navigate to the folder. (Also note
that Linux uses ``/`` instead of ``\`` for file paths.)

This method provides better compatibility than Cygwin or MSYS2 and more
convenience than a virtual machine if you partially need a Linux environment
for your workflow, however this method is only viable for Windows 10 64-bit
users.

2. Build the Documentation
--------------------------

The documentation can be built by running the ``makefile`` in the ``doc``
subdirectory.

To start, in your preferred web browser, use the drop down menu and select
“open file” to navigate into the sympy/doc folder saved on your computer. In
the doc folder, select the _build folder, then the html folder, and in the html
folder, open the index.html file.

To build the HTML documentation, run::

   cd doc

   make html

This builds a local version of the documentation in ``doc/_build/html`` in your
web browser.

Open ``_build/html/index.html``.

3. Make a Contribution
----------------------

For in-depth instructions on how to contribute to SymPy’s code base including
coding conventions, creating your environment, picking an issue to fix, and
opening a pull request, please read our full `Development Workflow
<https://github.com/sympy/sympy/wiki/Development-workflow>`_ guide.

Narrative Documentation Guidelines
==================================

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

.. _style_guide_docstring_guidelines:

Docstring Guidelines
====================

To contribute to SymPy’s docstrings, please read these guidelines in full.

A documentation string (docstring) is a string literal that occurs as the first
statement in a module, function, class, or method definition. Such a docstring
becomes the ``__doc__`` special attribute of that object.

**Example**

Here is a basic docstring::

    def fdiff(self, argindex=1):
        """
        Returns the first derivative of a Heaviside Function.

        Examples
        ========

        >>> from sympy import Heaviside, diff
        >>> from sympy.abc import x

        >>> Heaviside(x).fdiff()
        DiracDelta(x)

        >>> Heaviside(x**2 - 1).fdiff()
        DiracDelta(x**2 - 1)

        >>> diff(Heaviside(x)).fdiff()
        DiracDelta(x, 1)

        """

Every public function, class, method, and module should have a docstring that
describes what it does. Documentation that is specific to a function or class
in the module should be in the docstring of that function or class. The module
level docstring should discuss the purpose and scope of the module, and give a
high-level example of how to use the functions or classes in the module. A
module docstring is the docstring at the very top of the file, for example, the
docstring for `solvers.ode
<https://github.com/sympy/sympy/blob/85e684f782c71d247b13af71f2f134a9d894507e/sympy/solvers/ode.py>`_.

A public function is one that is intended to be used by the end-
user, or the public. Documentation is important for public functions because
they will be seen and used by many people.

A private function, on the other hand, is one that is only intended to be used
in the code in SymPy itself. Although it is less important to document private
functions, it also helps to have docstrings on private functions to help other
SymPy developers understand how to use the function.

It may not always be clear what is a public function and what is a private
function. If a function begins with an underscore, it is private, and if a
function is included in ``__init__.py`` it is public, but the converse is not
always true, so sometimes you have to decide based on context. In general, if
you are unsure, having documentation on a function is better than not having
documentation, regardless if it is public or private.

Docstrings should contain information aimed at users of the function. Comments
specific to the code or other notes that would only distract users should go in
comments in the code, not in docstrings.

Every docstring should have examples that show how the function works. Examples
are the most important part of a docstring. A single example showing input and
output to a function can be more helpful than a paragraph of descriptive text.

Remember that the primary consumers of docstrings are other human beings, not
machines, so it is important to describe what the function does in plain
English. Likewise, examples of how to use the function should be designed for
human readers, not just for the doctest machinery.

Keep in mind that while Sphinx is the primary way users consume docstrings, and
therefore the first platform to keep in mind while writing docstrings
(especially for public functions), it is not the only way users consume
docstrings. You can also view docstrings using ``help()`` or ``?`` in IPython.
When using ``help()``, for instance, it will show you all of the docstrings on
private methods. Additionally, anyone reading the source code directly will see
every docstring.

All public functions, classes, and methods and their corresponding docstrings
should be imported into the Sphinx docs, instructions on which can be found at
the end of this guide.

.. _style_guide_docstring_formatting:

Docstring Formatting
--------------------

Docstrings are are written in `reStructuredText
<http://docutils.sourceforge.net/rst.html>`_ format extended by `Sphinx
<http://www.sphinx-doc.org/en/master/>`_. Here is a concise guide to `Quick
reStructuredText <http://docutils.sourceforge.net/docs/user/rst/quickref.html>`_. More in-depth
information about using reStructuredText can be found in the `Sphinx
Documentation
<http://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html>`_.

In order for Sphinx to render docstrings nicely in the HTML documentation, some
formatting guidelines should be followed when writing docstrings:

* Always use """triple double quotes""" around docstrings. Use r"""raw triple
  double quotes""" if you use any backslashes in your docstrings.
* Include a blank line before the docstring’s closing quotes.
* Lines should not be longer than 80 characters.
* Always write class-level docstrings under the class definition line, as that
  is better to read in the source code.
* The various methods on the class can be mentioned in the docstring or
  examples if they are important, but details on them should go in the
  docstring for the method itself.
* Be aware that :: creates code blocks, which are rarely used in the
  docstrings. Any code example with example Python should be put in a doctest.
  Always check that the final version as rendered by Sphinx looks correct in
  the HTML.
* In order to make section underlining work nicely in docstrings, `numpydoc
  Sphinx extension <https://pypi.org/project/numpydoc/>`_ is used.
* Always double check that you have formatted your docstring correctly:

1. Make sure that your docstring is imported into Sphinx.
2. Build the Sphinx docs (``cd doc; make html``).
3. Make sure that Sphinx doesn't output any errors.
4. Open the page in ``_build/html`` and make sure that it is formatted
   correctly.

Docstring Sections
------------------

In SymPy’s docstrings, it is preferred that function, class, and method
docstrings consist of the following sections in this order:

1. Single-Sentence Summary
2. Explanation
3. Examples
4. Parameters
5. See Also
6. References

The Single-Sentence Summary and Examples sections are **required** for every
docstring. Docstrings will not pass review if these sections are not included.

Do not change the names of these supported sections, for example, the heading
“Examples” as a plural should be used even if there is only one example.

SymPy will continue to support all of the section headings listed in the `NumPy
Docstring Guide <https://numpydoc.readthedocs.io/en/latest/format.html>`_.

Headings should be underlined with the same length in equals signs.

If a section is not required and that information for the function in question
is unnecessary, do not use it. Unnecessary sections and cluttered docstrings
can make a function harder to understand. Aim for the minimal amount of
information required to understand the function.

1. Single-Sentence Summary
--------------------------

This section is **required** for every docstring. A docstring will not pass
review if it is not included. No heading is necessary for this section.

This section consists of a one-line sentence ending in a period that describes
the function, class, or method's effect.

Deprecation warnings should go directly after the Single-Sentence Summary, so
as to notify users right away. Deprecation warnings should be written as a note
in the Sphinx directive::

    .. note:: Deprecated in Sympy 0.7.1.

2. Explanation Section
----------------------

This section is encouraged. If you choose to include an Explanation section in
your docstring, it should be labeled with the heading “Explanation” underlined
with the same length in equals signs.

::

    Explanation
    ===========

This section consists of a more elaborate description of what the function,
class, or method does when the concise Single-Sentence Summary will not
suffice. This section should be used to clarify functionality in several
sentences or paragraphs.

3. Examples Section
-------------------

This section is **required** for every docstring. A docstring will not pass
review if it is not included. It should be labeled with the heading “Examples”
(even if there is only one example) underlined with the same length in equals
signs.

::

    Examples
    ========

This section consists of examples that show how the function works, called
doctests. Doctests should be complicated enough that they fully demonstrate the
API and functionality of the function, but simple enough that a user can
understand them without too much thought. The perfect doctest tells the user
exactly what they need to know about the function without reading any other
part of the docstring.

There should always be a blank line before the doctest. When multiple examples
are provided, they should be separated by blank lines. Comments explaining the
examples should have blank lines both above and below them.

Do not think of doctests as tests. Think of them as examples that happen to be
tested. They should demonstrate the API of the function to the user (i.e., what
the input parameters look like, what the output looks like, and what it does).
If you only want to test something, add a test to the relevant ``test_*.py file``.

You can use the ``./bin/coverage_doctest.py`` script to test the doctest
coverage of a file or module. Run the doctests with ``./bin/doctest``.

You should only skip the testing of an example if it is impossible to test it.
If necessary, testing of an example can be skipped by adding a special comment.

**Example**

.. code::

    >>> import random
    >>> random.random()      # doctest: +SKIP
    0.6868680200532414

If an example is longer than 80 characters, it should be line wrapped. Examples
should be line wrapped so that they are still valid Python code, using ``...``
continuation as in a Python prompt. For example, from the ODE module
documentation:

**Example**

.. code::

    >>> from sympy import Function, dsolve, cos, sin
    >>> from sympy.abc import x
    >>> f = Function('f')
    >>> dsolve(cos(f(x)) - (x*sin(f(x)) - f(x)**2)*f(x).diff(x),
    ... f(x), hint='1st_exact')
    Eq(x*cos(f(x)) + f(x)**3/3, C1)

Here ``dsolve(cos(f(x)) - (x*sin(f(x)) - f(x)**2)*f(x).diff(x), f(x), hint='1st_exact')`` is too long, so we line break it after a comma so that it
is readable, and put ``...`` on the continuation lines. If this is not done
correctly, the doctests will fail.

The output of a command can also be line wrapped. No ``...`` should be used in
this case. The doctester automatically accepts output that is line wrapped.

**Example**

.. code::

    >>> list(range(30))
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
    21, 22, 23, 24, 25, 26, 27, 28, 29]

In a doctest, write imports like ``sympy import ...`` instead of ``import
sympy`` or ``from sympy import *``. To define symbols, use ``from sympy.abc
import x``, unless the name is not in ``sympy.abc`` (for instance, if it has
assumptions), in which case use ``symbols`` like ``x, y = symbols('x y')``.

In general, you should run ``./bin/doctest`` to make sure your examples run
correctly, and fix them if they do not.

4. Parameters Section
---------------------

This section is encouraged. If you choose to include a Parameters section in
your docstring, it should be labeled with the heading “Parameters” underlined
with the same length in equals signs.

::

    Parameters
    ==========

If you have parameters listed in parentheses after a function, class, or method
name, you must include a parameters section.

This section consists of descriptions of the function arguments, keywords, and
their respective types.

Enclose variables in double backticks. The colon must be preceded by a space,
or omitted if the type is absent. For the parameter types, be as precise as
possible. If it is not necessary to specify a keyword argument, use
``optional``. Optional keyword parameters have default values, which are
displayed as part of the function signature. They can also be detailed in the
description.

When a parameter can only assume one of a fixed set of values, those values can
be listed in braces, with the default appearing first. When two or more input
parameters have exactly the same type, shape, and description, they can be
combined.

If the Parameters section is not formatted correctly, the documentation build
will render incorrectly.

If you wish to include a Returns section, write it as its own section with its
own heading.

**Example**

Here is an example of a correctly formatted Parameters section::

    def opt_cse(exprs, order='canonical'):
        """
        Find optimization opportunities in Adds, Muls, Pows and negative
        coefficient Muls.

        Parameters
        ==========

        exprs : list of sympy expressions
            The expressions to optimize.
        order : string, 'none' or 'canonical'
            The order by which Mul and Add arguments are processed. For large
            expressions where speed is a concern, use the setting order='none'.

        """

.. _style_guide_see_also:

5. See Also Section
-------------------

This section is encouraged. If you choose to include a See Also section in your
docstring, it should be labeled with the heading “See Also” underlined with the
same length in equals signs.

::

    See Also
    ========

This section consists of a listing of related functions, classes, and methods.
The related items can be described with a concise fragment (not a full
sentence) if desired, but this is not required. If the description spans more
than one line, subsequent lines must be indented.

The See Also section should only be used to reference other SymPy objects.
Anything that is a link should be embedded as a hyperlink in the text of the
docstring instead; see the References section for details.

Do not reference classes with ``class:Classname``, ``class:`Classname```, or
``:class:`Classname```, but rather only by their class name.

**Examples**

Here is a correctly formatted See Also section with concise descriptions::

    class erf(Function):
        r"""
        The Gauss error function.

        See Also
        ========

        erfc: Complementary error function.
        erfi: Imaginary error function.
        erf2: Two-argument error function.
        erfinv: Inverse error function.
        erfcinv: Inverse Complementary error function.
        erf2inv: Inverse two-argument error function.

        """

Here is a correctly formatted See Also section with just a list of names::

    class besselj(BesselBase):
        r"""
        Bessel function of the first kind.

        See Also
        ========

        bessely, besseli, besselk

        """

6. References Section
---------------------

This section is encouraged. If you choose to include a References section in
your docstring, it should be labeled with the heading “References” underlined
with the same length in equals signs.

::

    References
    ==========

This section consists of a list of references cited anywhere in the previous
sections. Any reference to other SymPy objects should go in the See Also
section instead.

The References section should include online resources, paper citations, and/or
any other printed resource giving general information about the function.
References are meant to augment the docstring, but should not be required to
understand it. References are numbered, starting from one, in the order in
which they are cited.

For online resources, only link to freely accessible and stable online
resources such as Wikipedia, Wolfram MathWorld, and the NIST Digital Library of
Mathematical Functions (DLMF), which are unlikely to suffer from hyperlink rot.

References for papers should include, in this order: reference citation, author
name, title of work, journal or publication, year published, page number.

If there is a DOI (Digital Object Identifier), include it in the citation and
make sure it is a clickable hyperlink.

**Examples**

Here is a References section that cites a printed resource::

    References
    ==========

    .. [1] [Kozen89] D. Kozen, S. Landau, Polynomial Decomposition Algorithms,
           Journal of Symbolic Computation 7 (1989), pp. 445-456

Here is a References section that cites printed and online resources::

    References
    ==========

    .. [1] Abramowitz, Milton; Stegun, Irene A., "Chapter 9," Handbook of
           Mathematical Functions with Formulas, Graphs, and Mathematical
           Tables, eds. (1965)
    .. [2] Luke, Y. L., The Special Functions and Their Approximations,
           Volume 1, (1969)
    .. [3] https://en.wikipedia.org/wiki/Bessel_function
    .. [4] http://functions.wolfram.com/Bessel-TypeFunctions/BesselJ/

Sample Docstring
================

Here is an example of a correctly formatted docstring::

    class gamma(Function):
        r"""
        The gamma function

        .. math::
           \Gamma(x) := \int^{\infty}_{0} t^{x-1} e^{-t} \mathrm{d}t.

        Explanation
        ===========

        The ``gamma`` function implements the function which passes through the
        values of the factorial function (i.e., $\Gamma(n) = (n - 1)!$), when n
        is an integer. More generally, $\Gamma(z)$ is defined in the whole
        complex plane except at the negative integers where there are simple
        poles.

        Examples
        ========

        >>> from sympy import S, I, pi, oo, gamma
        >>> from sympy.abc import x

        Several special values are known:

        >>> gamma(1)
        1
        >>> gamma(4)
        6
        >>> gamma(S(3)/2)
        sqrt(pi)/2

        The ``gamma`` function obeys the mirror symmetry:

        >>> from sympy import conjugate
        >>> conjugate(gamma(x))
        gamma(conjugate(x))

        Differentiation with respect to $x$ is supported:

        >>> from sympy import diff
        >>> diff(gamma(x), x)
        gamma(x)*polygamma(0, x)

        Series expansion is also supported:

        >>> from sympy import series
        >>> series(gamma(x), x, 0, 3)
        1/x - EulerGamma + x*(EulerGamma**2/2 + pi**2/12) + x**2*(-EulerGamma*pi**2/12 +
        polygamma(2, 1)/6 - EulerGamma**3/6) + O(x**3)

        We can numerically evaluate the ``gamma`` function to arbitrary
        precision on the whole complex plane:

        >>> gamma(pi).evalf(40)
        2.288037795340032417959588909060233922890
        >>> gamma(1+I).evalf(20)
        0.49801566811835604271 - 0.15494982830181068512*I

        See Also
        ========

        lowergamma: Lower incomplete gamma function.
        uppergamma: Upper incomplete gamma function.
        polygamma: Polygamma function.
        loggamma: Log Gamma function.
        digamma: Digamma function.
        trigamma: Trigamma function.
        beta: Euler Beta function.

        References
        ==========

        .. [1] https://en.wikipedia.org/wiki/Gamma_function
        .. [2] http://dlmf.nist.gov/5
        .. [3] http://mathworld.wolfram.com/GammaFunction.html
        .. [4] http://functions.wolfram.com/GammaBetaErf/Gamma/

        """

Docstrings for Classes that are Mathematical Functions
======================================================

SymPy is unusual in that it also has classes that are mathematical functions.
The docstrings for classes that are mathematical functions should include
details specific to this kind of class, as noted below:

* The Explanation section should include a mathematical definition of the
  function. This should use LaTeX math. Use $$ for :ref:`inline math
  <style_guide_math_formatting>` and .. math:: for display math, which should be
  used for the main definition. The variable names in the formulas should match
  the names of the parameters, and the LaTeX formatting should match the LaTeX
  pretty printing used by SymPy. As relevant, the mathematical definitions
  should mention their domain of definition, especially if the domain is
  different from the complex numbers.

* If there are multiple conventions in the literature for a function, make sure
  to clearly specify which convention SymPy uses.

* The Explanation section may also include some important mathematical facts
  about the function. These can alternately be demonstrated in the Examples
  section. Mathematical discussions should not be too long, as users can check
  the references for more details.

* The docstring does not need to discuss every implementation detail such as at
  which operations are defined on the function or at which points it evaluates
  in the "eval" method. Important or illuminating instances of these can be
  shown in the Examples section.

* The docstring should go on the class level (right under the line that has
  "class"). The "eval" method should not have a docstring.

* Private methods on the class, that is, any method that starts with an
  underscore, do not need to be documented. They can still be documented if you
  like, but note that these docstrings are not pulled into the Sphinx
  documentation, so they will only be seen by developers who are reading the
  code, so if there is anything very important that you want to mention here,
  it should go in the class-level docstring as well.

Best Practices for Writing Docstrings
=====================================

When writing docstrings, please follow all of the same formatting, style, and
tone preferences as when writing narrative documentation. For guidelines, see
:ref:`Best Practices for Writing Documentation
<style_guide_best_practices_for_writing_documentation>`, Formatting, Style, and
Tone.

Importing Docstrings into the Sphinx Documentation
==================================================

Here is a part of the ``doc/src/modules/geometry.txt`` file that imports the
relevant docstrings from geometry module into documentation::

    Utils
    =====

    .. module:: sympy.geometry.util

    .. autofunction:: intersection

    .. autofunction:: convex_hull

    .. autofunction:: are_similar

    Points
    ======

    .. module:: sympy.geometry.point

    .. autoclass:: Point
       :members:

    Lines
    =====

    .. module:: sympy.geometry.line

    .. autoclass:: LinearEntity
       :members:

    .. autoclass:: Line
       :members:

    .. autoclass:: Ray
       :members:

    .. autoclass:: Segment
       :members:

    Curves
    ======

    .. module:: sympy.geometry.curve

    .. autoclass:: Curve
       :members:

    Ellipses
    ========

    .. module:: sympy.geometry.ellipse

    .. autoclass:: Ellipse
       :members:

    .. autoclass:: Circle
       :members:

    Polygons
    ========

    .. module:: sympy.geometry.polygon

    .. autoclass:: Polygon
      :members:

    .. autoclass:: RegularPolygon
       :members:

    .. autoclass:: Triangle
       :members:

First namespace is set to particular submodule (file) with ``.. module::``
directive, then docstrings are imported with ``.. autoclass::`` or ``..
autofunction::`` relative to that submodule (file). Other methods are either
cumbersome to use (using full paths for all objects) or break something
(importing relative to main module using ``.. module:: sympy.geometry`` breaks
viewcode Sphinx extension). All files in ``doc/src/modules/`` should use this
format.

.. _style_guide_cross-referencing:

Cross-Referencing
=================

Any text that references another SymPy function should be formatted so that a
cross-reference link to that function's documentation is created automatically.
This is done using the RST cross-reference syntax. There are two different kinds
of objects that have conventions here:

1. Objects that are included in ``from sympy import *``, for example,
``sympy.acos``.

For these, use ``:obj:`~.acos()```. The ``~`` makes it so that
the text in the rendered HTML only shows ``acos``. Without it, it would use the
fully qualified name ``sympy.functions.elementary.trigonometric.acos``. However,
for names that are part of the global ``sympy`` namespace, we do not want to
encourage accessing them from their specific submodule, as this is an
implementation detail that could change. The ``.`` makes it so that the function
name is found automatically. Sometimes, Sphinx will give a warning that there
are multiple names found. If that happens, replace the ``.`` with the full name.
For example, ``:obj:`~sympy.solvers.solvers.solve()```. For functions, methods,
and classes, it is a convention to add () after the name to indicate such.

You may also use a more specific type indicator instead of ``obj`` (see
https://www.sphinx-doc.org/en/master/usage/restructuredtext/domains.html#cross-referencing-python-objects).
However, ``obj`` will always work, and sometimes SymPy names are not the type
you might expect them to be. For example, mathematical function objects such as
``sin`` are not actually a Python function, rather they are a Python class,
therefore ``:func:`~.sin``` will not work.

2. Objects that are not included in ``from sympy import *``, for example,
``sympy.physics.vector.dynamicsymbols``.

This can be public API objects from submodules that are not included in the main
``sympy/__init__.py``, such as the physics submodule, or private API objects
that are not necessarily intended to be used by end-users (but should still be
documented). In this case, you must show the fully qualified name, so do not use
the ``~.`` syntax. For example,
``:obj:`sympy.physics.vector.dynamicsymbols()```.

You may also write custom text that links to the documentation for something
using the following syntax ``:obj:`custom text<object>```. For example,
``:obj:`the sine function <.sin>``` produces the text "the sine function" that
links to the documentation for ``sin``. Note that the ``~`` character should
not be used here.

Note that references in the :ref:`See Also <style_guide_see_also>` section of
the docstrings do not require the ``:obj:`` syntax.

If the resulting cross reference is written incorrectly, Sphinx will error when
building the docs with an error like:

::

   WARNING: py:obj reference target not found: expand

Here are some troubleshooting tips to fix the errors:

* Make sure you have used the correct syntax, as described above.
* Make sure you spelled the function name correctly.
* Check if the function you are trying to cross-reference is actually included
  in the Sphinx documentation. If it is not, Sphinx will not be able to create
  a reference for it. In that case, you should add it to the appropriate RST
  file as described in the :ref:`Docstring Guidelines
  <style_guide_docstring_guidelines>`.
* If the function or object is not included in ``from sympy import
  *``, you will need to use the fully qualified name, like
  ``sympy.submodule.submodule.function`` instead of just ``function``.
* A fully qualified name must include the full submodule for a function all the
  way down to the file. For example, ``sympy.physics.vector.ReferenceFrame``
  will not work (even though you can access it that way in code). It has to be
  ``sympy.physics.vector.frame.ReferenceFrame``.
* If the thing you are referring to does not actually have somewhere to link
  to, do not use the ``:obj:`` syntax. Instead, mark it as code using double
  backticks. Examples of things that cannot be linked to are Python built in
  functions like ``int`` or ``NotImplementedError``, functions from other
  modules outside of SymPy like ``matplotlib.plot``, and variable or parameter
  names that are specific to the text at hand. In general, if the object cannot
  be accessed as ``sympy.something.something.object``, it cannot be cross-
  referenced and you should not use the ``:obj:`` syntax.
* If you are using are using one of the `type specific
  <https://www.sphinx-doc.org/en/master/usage/restructuredtext/domains.html#cross-referencing-python-objects>`_
  identifiers like ``:func:``, be sure that the type for it is correct.
  ``:func:`` only refers to Python functions. For classes, you need to use
  ``:class:``, and for methods on a class you need to use ``:method:``. In
  general, it is recommended to use ``:obj:``, as this will work for any type
  of object.
* If you cannot get the cross-referencing syntax to work, go ahead and submit
  the pull request as is and ask the reviewers for help.

You may also see errors like:

::

    WARNING: more than one target found for cross-reference 'subs()':
    sympy.core.basic.Basic.subs, sympy.matrices.common.MatrixCommon.subs,
    sympy.physics.vector.vector.Vector.subs,
    sympy.physics.vector.dyadic.Dyadic.subs

for instance, from using ``:obj:`~.subs```. This means that the ``.`` is not
sufficient to find the function, because there are multiple names in SymPy
named ``subs``. In this case, you need to use the fully qualified name. You can
still use ``~`` to make it shortened in the final text, like
``:obj:`~sympy.core.basic.Basic.subs```.

The line numbers for warnings in Python files are relative to the top of the
docstring, not the file itself. The line numbers are often not completely
correct, so you will generally have to search the docstring for the part that
the warning is referring to. This is due to a bug in Sphinx.
