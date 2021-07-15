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
