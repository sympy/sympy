For instructions on installation, building the documentation, and guidelines for
contributing to SymPy's documentation, please read the `SymPy Documentation
Style Guide <https://docs.sympy.org/dev/documentation-style-guide.html>`_.

The SymPy Documentation Style Guide can also be read at
src/documentation-style-guide.rst.


Speeding up local Sphinx builds
-------------------------------

While creating the documentation locally, Sphinx might take some time because of the
``sphinx.ext.viewcode`` extension, which creates links to the source code.

If you are working on your documentation and want the build process to be faster, you
can temporarily comment out ``'sphinx.ext.viewcode'`` from ``doc/src/conf.py``.

This is meant for local development and should not be committed.
