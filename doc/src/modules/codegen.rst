.. _codegen_prose:

================================================
Structural Details of Code Generation with Sympy
================================================

Several submodules in SymPy allow one to generate directly compilable and
executable code in a variety of different programming languages from Sympy
expressions. In addition, there are functions that generate Python importable
objects that can evaluate SymPy expressions very efficiently.

We will start with a brief introduction to the components that make up the code
generation capabilities of SymPy.

Introduction
------------

There are four main levels of abstractions::

   expression
      |
   code printers
      |
   code generators
      |
   autowrap

:mod:`sympy.utilities.autowrap` uses codegen, and codegen uses the code
printers. :mod:`sympy.utilities.autowrap` does everything: it lets you go
from SymPy expression to numerical function in the same Python process in one
step. codegen is actual code generation, i.e., to compile and use later, or to
include in some larger project.

The code printers translate the SymPy objects into actual code, like ``abs(x)
-> fabs(x)`` (for C).

The code printers don't print optimal code in many cases. An example of this is
powers in C. ``x**2`` prints as ``pow(x, 2)`` instead of ``x*x``.  Other
optimizations (like mathematical simplifications) should happen before the code
printers.

Currently, :py:func:`sympy.simplify.cse_main.cse` is not applied automatically anywhere in this
chain. It ideally happens at the codegen level, or somewhere above it.

We will iterate through the levels below.

The following three lines will be used to setup each example::

    >>> from sympy import *
    >>> init_printing(use_unicode=True)
    >>> from sympy.abc import a, e, k, n, r, t, x, y, z, T, Z
    >>> from sympy.abc import beta, omega, tau
    >>> f, g = symbols('f, g', cls=Function)

Code printers (sympy.printing)
------------------------------

This is where the meat of code generation is; the translation of SymPy
expressions to specific languages. Supported languages are C
(:py:func:`sympy.printing.ccode.ccode`), Fortran 95
(:py:func:`sympy.printing.fcode.fcode`), Javascript
(:py:func:`sympy.printing.jscode.jscode`), Julia
(:py:func:`sympy.printing.julia.julia_code`), Mathematica
(:py:func:`sympy.printing.mathematica.mathematica_code`), Octave/Matlab
(:py:func:`sympy.printing.octave.octave_code`), Python (print_python, which is
actually more like a lightweight version of codegen for Python, and
:py:func:`sympy.printing.lambdarepr.lambdarepr`, which supports many libraries
(like NumPy), and theano
(:py:func:`sympy.printing.theanocode.theano_function`). The code printers are
special cases of the other prints in SymPy (str printer, pretty printer, etc.).

An important distinction is that the code printer has to deal with assignments
(using the :class:`sympy.printing.codeprinter.Assignment` object).This serves as
building blocks for the code printers and hence the ``codegen`` module.  An
example that shows the use of ``Assignment``::

    >>> from sympy.printing.codeprinter import Assignment
    >>> mat = Matrix([x, y, z]).T
    >>> known_mat = MatrixSymbol('K', 1, 3)
    >>> Assignment(known_mat, mat)
    K := [x  y  z]
    >>> Assignment(known_mat, mat).lhs
    K
    >>> Assignment(known_mat, mat).rhs
    [x  y  z]

Here is a simple example of printing a C version of a SymPy expression::

    >>> expr = (Rational(-1, 2) * Z * k * (e**2) / r)
    >>> expr
        2
    -Z⋅e ⋅k
    ────────
      2⋅r
    >>> ccode(expr)
    -1.0L/2.0L*Z*pow(e, 2)*k/r
    >>> ccode(expr, assign_to="E")
    E = -1.0L/2.0L*Z*pow(e, 2)*k/r;

``Piecewise`` expressions are converted into conditionals. If an ``assign_to``
variable is provided an if statement is created, otherwise the ternary operator
is used. Note that if the ``Piecewise`` lacks a default term, represented by
``(expr, True)`` then an error will be thrown.  This is to prevent generating
an expression that may not evaluate to anything. A use case for ``Piecewise``::

    >>> expr = Piecewise((x + 1, x > 0), (x, True))
    >>> print(fcode(expr, tau))
          if (x > 0) then
             tau = x + 1
          else
             tau = x
          end if

The various printers also tend to support ``Indexed`` objects well. With
``contract=True`` these expressions will be turned into loops, whereas
``contract=False`` will just print the assignment expression that should be
looped over::

    >>> len_y = 5
    >>> mat_1 = IndexedBase('mat_1', shape=(len_y,))
    >>> mat_2 = IndexedBase('mat_2', shape=(len_y,))
    >>> Dy = IndexedBase('Dy', shape=(len_y-1,))
    >>> i = Idx('i', len_y-1)
    >>> eq = Eq(Dy[i], (mat_1[i+1] - mat_1[i]) / (mat_2[i+1] - mat_2[i]))
    >>> print(jscode(eq.rhs, assign_to=eq.lhs, contract=False))
    Dy[i] = (mat_1[i + 1] - mat_1[i])/(mat_2[i + 1] - mat_2[i]);
    >>> Res = IndexedBase('Res', shape=(len_y,))
    >>> j = Idx('j', len_y)
    >>> eq = Eq(Res[j], mat_1[j]*mat_2[j])
    >>> print(jscode(eq.rhs, assign_to=eq.lhs, contract=True))
    for (var j=0; j<5; j++){
       Res[j] = 0;
    }
    for (var j=0; j<5; j++){
       for (var j=0; j<5; j++){
          Res[j] = Res[j] + mat_1[j]*mat_2[j];
       }
    }
    >>> print(jscode(eq.rhs, assign_to=eq.lhs, contract=False))
    Res[j] = mat_1[j]*mat_2[j];


Custom printing can be defined for certain types by passing a dictionary of
"type" : "function" to the ``user_functions`` kwarg. Alternatively, the
dictionary value can be a list of tuples i.e., ``[(argument_test,
cfunction_string)]``. This can be used to call a custom Octave function::

    >>> custom_functions = {
    ...   "f": "existing_octave_fcn",
    ...   "g": [(lambda x: x.is_Matrix, "my_mat_fcn"),
    ...         (lambda x: not x.is_Matrix, "my_fcn")]
    ... }
    >>> mat = Matrix([[1, x]])
    >>> octave_code(f(x) + g(x) + g(mat), user_functions=custom_functions)
    existing_octave_fcn(x) + my_fcn(x) + my_mat_fcn([1 x])

An example of Mathematica code printer::

    >>> x_ = Function('x')
    >>> expr = x_(n*T) * sin((t - n*T) / T)
    >>> expr = expr / ((-T*n + t) / T)
    >>> expr
                ⎛-T⋅n + t⎞
    T⋅x(T⋅n)⋅sin⎜────────⎟
                ⎝   T    ⎠
    ──────────────────────
           -T⋅n + t

    >>> expr = summation(expr, (n, -1, 1))
    >>> mathematica_code(expr)
    T*x[-T]*Sin[(T + t)/T]/(T + t) + T*x[T]*Sin[(-T + t)/T]/(-T + t) + T*x[0]*Sin[
    t/T]/t

We can go through a common expression in different languages we support and see
how it works::

    >>> k, g1, g2, r, I, S = symbols("k, gamma_1, gamma_2, r, I, S")
    >>> expr = k * g1 * g2 / (r**3)
    >>> expr = expr * 2 * I * S * (3 * (cos(beta))**2 - 1) / 2
    >>> expr
                ⎛     2       ⎞
    I⋅S⋅γ₁⋅γ₂⋅k⋅⎝3⋅cos (β) - 1⎠
    ───────────────────────────
                  3
                 r
    >>> print(jscode(expr, assign_to="H_is"))
    H_is = I*S*gamma_1*gamma_2*k*(3*Math.pow(Math.cos(beta), 2) - 1)/Math.pow(r, 3);
    >>> print(ccode(expr, assign_to="H_is"))
    H_is = I*S*gamma_1*gamma_2*k*(3*pow(cos(beta), 2) - 1)/pow(r, 3);
    >>> print(fcode(expr, assign_to="H_is"))
          H_is = I*S*gamma_1*gamma_2*k*(3*cos(beta)**2 - 1)/r**3
    >>> print(julia_code(expr, assign_to="H_is"))
    H_is = I.*S.*gamma_1.*gamma_2.*k.*(3*cos(beta).^2 - 1)./r.^3
    >>> print(octave_code(expr, assign_to="H_is"))
    H_is = I.*S.*gamma_1.*gamma_2.*k.*(3*cos(beta).^2 - 1)./r.^3;
    >>> print(mathematica_code(expr))
    I*S*gamma_1*gamma_2*k*(3*Cos[beta]^2 - 1)/r^3

Codegen (sympy.utilities.codegen)
---------------------------------

This module deals with creating compilable code from SymPy expressions. This is
lower level than autowrap, as it doesn't actually attempt to compile the code,
but higher level than the printers, as it generates compilable files (including
header files), rather than just code snippets.

The user friendly functions, here, are ``codegen`` and ``make_routine``.
``codegen`` takes a list of ``(variable, expression)`` pairs and a language (C,
F95, and Octave/Matlab are supported). It returns, as strings, a code file and
a header file (for relevant languages). The variables are created as functions
that return the value of the expression as output.

.. note:: The ``codegen`` callable is not in the sympy namespace automatically,
   to use it you must first import ``codegen`` from ``sympy.utilities.codegen``

For instance::

    >>> from sympy.utilities.codegen import codegen
    >>> length, breadth, height = symbols('length, breadth, height')
    >>> [(c_name, c_code), (h_name, c_header)] = \
    ... codegen(('volume', length*breadth*height), "C", "test",
    ...         header=False, empty=False)
    >>> print(c_name)
    test.c
    >>> print(c_code)
    #include "test.h"
    #include <math.h>
    double volume(double breadth, double height, double length) {
       double volume_result;
       volume_result = breadth*height*length;
       return volume_result;
    }
    >>> print(h_name)
    test.h
    >>> print(c_header)
    #ifndef PROJECT__TEST__H
    #define PROJECT__TEST__H
    double volume(double breadth, double height, double length);
    #endif

Various flags to ``codegen`` let you modify things. The project name for
preprocessor instructions can be varied using ``project``. Variables listed as
global variables in arg ``global_vars`` will not show up as function arguments.

``language`` is a case-insensitive string that indicates the source code
language. Currently, ``C``, ``F95`` and ``Octave`` are supported. ``Octave``
generates code compatible with both Octave and Matlab.

``header`` when True, a header is written on top of each source file. ``empty``
when True, empty lines are used to structure the code. With
``argument_sequence`` a sequence of arguments for the routine can be defined in
a preferred order.

``prefix`` defines a prefix for the names of the files that contain the source
code.  If omitted, the name of the first name_expr tuple is used.

``to_files`` when True, the code will be written to one or more files with the
given prefix.

Here is an example::

    >>> [(f_name, f_code), header] = codegen(("volume", length*breadth*height),
    ...     "F95", header=True, empty=False, argument_sequence=(breadth, length),
    ...     global_vars=(height,))
    >>> print(f_code)
    !******************************************************************************
    !*                    Code generated with sympy ...                           *
    !*                                                                            *
    !*              See http://www.sympy.org/ for more information.               *
    !*                                                                            *
    !*                       This file is part of 'project'                       *
    !******************************************************************************
    REAL*8 function volume(breadth, length)
    implicit none
    REAL*8, intent(in) :: breadth
    REAL*8, intent(in) :: length
    volume = breadth*height*length
    end function

The method ``make_routine`` creates a ``Routine`` object, which represents an
evaluation routine for a set of expressions. This is only good for internal use
by the CodeGen objects, as an intermediate representation from SymPy expression
to generated code.  It is not recommended to make a ``Routine`` object
yourself. You should instead use ``make_routine`` method. ``make_routine`` in
turn calls the ``routine`` method of the CodeGen object depending upon the
language of choice. This creates the internal objects representing assignments
and so on, and creates the ``Routine`` class with them.

The various codegen objects such as ``Routine`` and ``Variable`` aren't SymPy
objects (they don't subclass from Basic).

For example::

    >>> from sympy.utilities.codegen import make_routine
    >>> from sympy.physics.hydrogen import R_nl
    >>> expr = R_nl(3, y, x, 6)
    >>> routine = make_routine('my_routine', expr)
    >>> [arg.result_var for arg in routine.results]   # doctest: +SKIP
    [result₅₁₄₂₃₄₁₆₈₁₃₉₇₇₁₉₄₂₈]
    >>> [arg.expr for arg in routine.results]
    ⎡                ___________                                           ⎤
    ⎢          y    ╱ (-y + 2)!   -2⋅x                                     ⎥
    ⎢4⋅√6⋅(4⋅x) ⋅  ╱  ───────── ⋅ℯ    ⋅assoc_laguerre(-y + 2, 2⋅y + 1, 4⋅x)⎥
    ⎢            ╲╱    (y + 3)!                                            ⎥
    ⎢──────────────────────────────────────────────────────────────────────⎥
    ⎣                                  3                                   ⎦
    >>> [arg.name for arg in routine.arguments]
    [x, y]

Another more complicated example with a mixture of specified and
automatically-assigned names.  Also has Matrix output::

    >>> routine = make_routine('fcn', [x*y, Eq(a, 1), Eq(r, x + r), Matrix([[x, 2]])])
    >>> [arg.result_var for arg in routine.results]   # doctest: +SKIP
    [result_5397460570204848505]
    >>> [arg.expr for arg in routine.results]
    [x⋅y]
    >>> [arg.name for arg in routine.arguments]   # doctest: +SKIP
    [x, y, a, r, out_8598435338387848786]

We can examine the various arguments more closely::

    >>> from sympy.utilities.codegen import (InputArgument, OutputArgument,
    ...                                      InOutArgument)
    >>> [a.name for a in routine.arguments if isinstance(a, InputArgument)]
    [x, y]

    >>> [a.name for a in routine.arguments if isinstance(a, OutputArgument)]  # doctest: +SKIP
    [a, out_8598435338387848786]
    >>> [a.expr for a in routine.arguments if isinstance(a, OutputArgument)]
    [1, [x  2]]

    >>> [a.name for a in routine.arguments if isinstance(a, InOutArgument)]
    [r]
    >>> [a.expr for a in routine.arguments if isinstance(a, InOutArgument)]
    [r + x]

The full API reference can be viewed :ref:`here<codegen_API>`.

Autowrap
--------

Autowrap automatically generates code, writes it to disk, compiles it, and
imports it into the current session. Main functions of this module are
``autowrap``, ``binary_function``, and ``ufuncify``.

It also automatically converts expressions containing ``Indexed`` objects into
summations. The classes IndexedBase, Indexed and Idx represent a matrix element
M[i, j]. See :ref:`tensor_module` for more on this.

.. _autowrap:

``autowrap`` creates a wrapper using f2py or Cython and creates a numerical function.

.. note:: The ``autowrap`` callable is not in the sympy namespace automatically,
   to use it you must first import ``autowrap`` from ``sympy.utilities.autowrap``

The callable returned from autowrap() is a binary Python function, not a SymPy
object. For example::

    >>> from sympy.utilities.autowrap import autowrap
    >>> expr = ((x - y + z)**(13)).expand()
    >>> binary_func = autowrap(expr)    # doctest: +SKIP
    >>> binary_func(1, 4, 2)    # doctest: +SKIP
    -1.0

The various flags available with autowrap() help to modify the services
provided by the method. The argument ``tempdir`` tells autowrap to compile the
code in a specific directory, and leave the files intact when finished. For
instance::

    >>> from sympy.utilities.autowrap import autowrap
    >>> from sympy.physics.qho_1d import psi_n
    >>> x_ = IndexedBase('x')
    >>> y_ = IndexedBase('y')
    >>> m = symbols('m', integer=True)
    >>> i = Idx('i', m)
    >>> qho = autowrap(Eq(y_[i], psi_n(0, x_[i], m, omega)), tempdir='/tmp')  # doctest: +SKIP

Checking the Fortran source code in the directory specified reveals this::

    subroutine autofunc(m, omega, x, y)
    implicit none
    INTEGER*4, intent(in) :: m
    REAL*8, intent(in) :: omega
    REAL*8, intent(in), dimension(1:m) :: x
    REAL*8, intent(out), dimension(1:m) :: y
    INTEGER*4 :: i

    REAL*8, parameter :: hbar = 1.05457162d-34
    REAL*8, parameter :: pi = 3.14159265358979d0
    do i = 1, m
       y(i) = (m*omega)**(1.0d0/4.0d0)*exp(-4.74126166983329d+33*m*omega*x(i &
             )**2)/(hbar**(1.0d0/4.0d0)*pi**(1.0d0/4.0d0))
    end do

    end subroutine

Using the argument ``args`` along with it changes argument sequence::

    >>> eq = Eq(y_[i], psi_n(0, x_[i], m, omega))
    >>> qho = autowrap(eq, tempdir='/tmp', args=[y, x, m, omega])  # doctest: +SKIP

yields::

    subroutine autofunc(y, x, m, omega)
    implicit none
    INTEGER*4, intent(in) :: m
    REAL*8, intent(in) :: omega
    REAL*8, intent(out), dimension(1:m) :: y
    REAL*8, intent(in), dimension(1:m) :: x
    INTEGER*4 :: i

    REAL*8, parameter :: hbar = 1.05457162d-34
    REAL*8, parameter :: pi = 3.14159265358979d0
    do i = 1, m
       y(i) = (m*omega)**(1.0d0/4.0d0)*exp(-4.74126166983329d+33*m*omega*x(i &
             )**2)/(hbar**(1.0d0/4.0d0)*pi**(1.0d0/4.0d0))
    end do

    end subroutine

The argument ``verbose`` is boolean, optional and if True, autowrap will not
mute the command line backends. This can be helpful for debugging.

The argument ``language`` and ``backend`` are used to change defaults:
``Fortran`` and ``f2py`` to ``C`` and ``Cython``. The argument helpers is used
to define auxiliary expressions needed for the main expression. If the main
expression needs to call a specialized function it should be put in the
``helpers`` iterable. Autowrap will then make sure that the compiled main
expression can link to the helper routine. Items should be tuples with
``(<function_name>, <sympy_expression>, <arguments>)``. It is mandatory to
supply an argument sequence to helper routines.

.. _binary_function:

Another method available at the ``autowrap`` level is ``binary_function``. It
returns a sympy function. The advantage is that we can have very fast functions
as compared to SymPy speeds. This is because we will be using compiled
functions with Sympy attributes and methods. An illustration::

    >>> from sympy.utilities.autowrap import binary_function
    >>> from sympy.physics.hydrogen import R_nl
    >>> psi_nl = R_nl(1, 0, a, r)
    >>> f = binary_function('f', psi_nl)    # doctest: +SKIP
    >>> f(a, r).evalf(3, subs={a: 1, r: 2})  # doctest: +SKIP
    0.766

.. _ufuncify_method:

While NumPy operations are very efficient for vectorized data but they
sometimes incur unnecessary costs when chained together.
Consider the following operation

    >>> x = get_numpy_array(...) # doctest: +SKIP
    >>> y = sin(x) / x

The operators ``sin`` and ``/`` call routines that execute tight for loops in
``C``. The resulting computation looks something like this

.. code:: c

    for(int i = 0; i < n; i++)
    {
        temp[i] = sin(x[i]);
    }
    for(int i = i; i < n; i++)
    {
        y[i] = temp[i] / x[i];
    }

This is slightly sub-optimal because

1.  We allocate an extra ``temp`` array
2.  We walk over ``x`` memory twice when once would have been sufficient

A better solution would fuse both element-wise operations into a single for loop

.. code:: c

    for(int i = i; i < n; i++)
    {
        y[i] = sin(x[i]) / x[i];
    }

Statically compiled projects like NumPy are unable to take advantage of such
optimizations. Fortunately, SymPy is able to generate efficient low-level C
or Fortran code. It can then depend on projects like ``Cython`` or ``f2py`` to
compile and reconnect that code back up to Python. Fortunately this process is
well automated and a SymPy user wishing to make use of this code generation
should call the ``ufuncify`` function.

``ufuncify`` is the third method available with Autowrap module. It basically
implies 'Universal functions' and follows an ideology set by Numpy. The main
point of ufuncify as compared to autowrap is that it allows arrays as arguments
and can operate in an element-by-element fashion. The core operation done
element-wise is in accordance to Numpy's array broadcasting rules. See `this
<http://docs.scipy.org/doc/numpy/reference/ufuncs.html>`_ for more.

    >>> from sympy import *
    >>> from sympy.abc import x
    >>> expr = sin(x)/x

    >>> from sympy.utilities.autowrap import ufuncify
    >>> f = ufuncify([x], expr) # doctest: +SKIP

This function ``f`` consumes and returns a NumPy array. Generally ``ufuncify``
performs at least as well as ``lambdify``. If the expression is complicated
then ``ufuncify`` often significantly outperforms the NumPy backed solution.
Jensen has a good `blog post <http://ojensen.wordpress.com/2010/08/10/fast-ufunc-ish-hydrogen-solutions/>`_ on this topic.

Let us see an example for some quantitative analysis::

    >>> from sympy.physics.hydrogen import R_nl
    >>> expr = R_nl(3, 1, x, 6)
    >>> expr
                    -2⋅x
    8⋅x⋅(-4⋅x + 4)⋅ℯ
    ────────────────────
             3

The lambdify function translates SymPy expressions into Python functions,
leveraging a variety of numerical libraries. By default lambdify relies on
implementations in the ``math`` standard library. Naturally, Raw Python is
faster than Sympy. However it also supports ``mpmath`` and most notably,
``numpy``. Using the numpy library gives the generated function access to
powerful vectorized ufuncs that are backed by compiled C code.

Let us compare the speeds::

    >>> from sympy.utilities.autowrap import ufuncify
    >>> from sympy.utilities.lambdify import lambdify
    >>> fn_numpy = lambdify(x, expr, 'numpy')   # doctest: +SKIP
    >>> fn_fortran = ufuncify([x], expr, backend='f2py')    # doctest: +SKIP
    >>> from numpy import linspace  # doctest: +SKIP
    >>> xx = linspace(0, 1, 5)  # doctest: +SKIP
    >>> fn_numpy(xx)    # doctest: +SKIP
    [ 0.          1.21306132  0.98101184  0.44626032  0.        ]
    >>> fn_fortran(xx)  # doctest: +SKIP
    [ 0.          1.21306132  0.98101184  0.44626032  0.        ]
    >>> import timeit
    >>> timeit.timeit('fn_numpy(xx)', 'from __main__ import fn_numpy, xx', number=10000)    # doctest: +SKIP
    0.18891601900395472
    >>> timeit.timeit('fn_fortran(xx)', 'from __main__ import fn_fortran, xx', number=10000)    # doctest: +SKIP
    0.004707066000264604

The options available with ufuncify are more or less the same as those
available with ``autowrap``.

There are other facilities available with Sympy to do efficient numeric
computation. See :ref:`this<numeric_computation>` page for a comparison among them.
