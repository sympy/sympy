===============
Autowrap Module
===============

The autowrap module works very well in tandem with the Indexed classes of the
:ref:`tensor_module`.  Here is a simple example that shows how to setup a binary
routine that calculates a matrix-vector product.

>>> from sympy.utilities.autowrap import autowrap
>>> from sympy import symbols, IndexedBase, Idx, Eq
>>> A, x, y = map(IndexedBase, ['A', 'x', 'y'])
>>> m, n = symbols('m n', integer=True)
>>> i = Idx('i', m)
>>> j = Idx('j', n)
>>> instruction = Eq(y[i], A[i, j]*x[j]); instruction
Eq(y[i], A[i, j]*x[j])

Because the code printers treat Indexed objects with repeated indices as a
summation, the above equality instance will be translated to low-level code for
a matrix vector product.  This is how you tell SymPy to generate the code,
compile it and wrap it as a python function:

>>> matvec = autowrap(instruction)                 # doctest: +SKIP

That's it.  Now let's test it with some numpy arrays.  The default wrapper
backend is f2py.  The wrapper function it provides is set up to accept python
lists, which it will silently convert to numpy arrays.  So we can test the
matrix vector product like this:

>>> M = [[0, 1],
...      [1, 0]]
>>> matvec(M, [2, 3])                              # doctest: +SKIP
[ 3.  2.]

Implementation details
======================

The autowrap module is implemented with a backend consisting of CodeWrapper
objects.  The base class ``CodeWrapper`` takes care of details about module
name, filenames and options.  It also contains the driver routine, which runs
through all steps in the correct order, and also takes care of setting up and
removing the temporary working directory.

The actual compilation and wrapping is done by external resources, such as the
system installed f2py command. The Cython backend runs a distutils setup script
in a subprocess. Subclasses of CodeWrapper takes care of these
backend-dependent details.

API Reference
=============

.. automodule:: sympy.utilities.autowrap
   :members:
