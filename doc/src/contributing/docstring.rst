.. _style_guide_docstring_guidelines:

=========================
Docstrings Style Guide
=========================

General Guidelines
--------------------

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

A public function is one that is intended to be used by the end-user,
or the public. Documentation is important for public functions because
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

Formatting
-------------

Docstrings are written in `reStructuredText
<https://docutils.sourceforge.io/rst.html>`_ format extended by `Sphinx
<https://www.sphinx-doc.org/en/master/>`_. Here is a concise guide to `Quick
reStructuredText <https://docutils.sourceforge.io/docs/user/rst/quickref.html>`_. More in-depth
information about using reStructuredText can be found in the `Sphinx
Documentation
<https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html>`_.

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

Sections
---------

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section is **required** for every docstring. A docstring will not pass
review if it is not included. No heading is necessary for this section.

This section consists of a one-line sentence ending in a period that describes
the function, class, or method's effect.

Deprecation warnings should go directly after the Single-Sentence Summary, so
as to notify users right away. Deprecation warnings should be written as a ``deprecated``
in the Sphinx directive::

    .. deprecated:: 1.1

       The ``simplify_this`` function is deprecated. Use :func:`simplify`
       instead. See its documentation for more information.

See :ref:`deprecation-documentation` for more details.

2. Explanation Section
^^^^^^^^^^^^^^^^^^^^^^^

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

.. _style_guide_docstring_examples_section:

3. Examples Section
^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^

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
    .. [4] https://functions.wolfram.com/Bessel-TypeFunctions/BesselJ/

Sample Docstring
------------------

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
        1/x - EulerGamma + x*(EulerGamma**2/2 + pi**2/12) +
        x**2*(-EulerGamma*pi**2/12 - zeta(3)/3 - EulerGamma**3/6) + O(x**3)

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
        .. [2] https://dlmf.nist.gov/5
        .. [3] https://mathworld.wolfram.com/GammaFunction.html
        .. [4] https://functions.wolfram.com/GammaBetaErf/Gamma/

        """

Docstrings for Classes that are Mathematical Functions
--------------------------------------------------------

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
---------------------------------------

When writing docstrings, please follow all of the same formatting, style, and
tone preferences as when writing narrative documentation. For guidelines, see
:ref:`Best Practices for Writing Documentation
<style_guide_best_practices_for_writing_documentation>`, Formatting, Style, and
Tone.

Importing Docstrings into the Sphinx Documentation
----------------------------------------------------

Here are excerpts from the ``doc/src/modules/geometry`` directory that imports the
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
------------------

Any text that references another SymPy function should be formatted so that a
cross-reference link to that function's documentation is created automatically.
This is done using the RST cross-reference syntax. There are two different kinds
of objects that have conventions here:

1. Objects that are included in ``from sympy import *``, for example,
``sympy.acos``.

For these, use ``:obj:`~.acos()```. The ``~`` makes it so that the text in the
rendered HTML only shows ``acos`` instead of the fully qualified name
``sympy.functions.elementary.trigonometric.acos``. (This will encourage importing
names from the global ``sympy`` namespace instead of a specific submodule.)
The ``.`` makes it so that the function name is found automatically. (If Sphinx gives
a warning that there are multiple names found, replace the ``.`` with
the full name.  For example, ``:obj:`~sympy.solvers.solvers.solve()```.) Adding a trailing
pair of parentheses is a convention for indicating the name is a function, method, or
class.

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
  be accessed as ``sympy.something.something.object``, it cannot be
  cross-referenced and you should not use the ``:obj:`` syntax.
* If you are using one of the `type specific
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
    sympy.core.basic.Basic.subs, sympy.matrices.matrixbase.MatrixBase.subs,
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
