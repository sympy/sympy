===============
Code Generation
===============
Several submodules in Sympy are provided to generate directly compilable 
code from Sympy expressions. Furthermore, there are submodules provided 
that help use other numeric systems as backend for construction and 
manipulation of mathematical expressions to speed up the process.

We will start with a brief introduction, and inspect the building 
blocks of the code generation and then the code generation itself.

Introduction
------------

There are three levels::

    autowrap
       |
    codegen
       |
    code printers

Autowrap uses codegen, and codegen uses the code printers. Autowrap does 
everything: it lets you go from SymPy expression to numerical 
function in the same Python process in one step. codegen is actual 
code generation, i.e., to compile and use later, or to include in some larger 
project.

The code printers translate the SymPy objects into actual code, 
like abs(x) -> fabs(x) (for C).

The code printers don't print optimal code in many cases. 
An example of this is powers in C. x**2 prints as pow(x, 2) instead of x*x. 
Other optimizations (like mathematical simplifications) should happen 
before the code printers.

cse is not applied anywhere in this chain. It ideally happens at the 
codegen level, or somewhere above it.

We will iterate through the levels, bottom up.

Code printers (sympy.printing)
------------------------------
This is where the meat of code generation is, the translation of SymPy
expressions to specific code.Supported languages are C (``ccode``), Fortran 95 (``fcode``), 
Javascript (``jscode``), Mathematica (``mathematica_code``), Octave/Matlab (``octave_code``), 
Python (print_python, which is actually more like a lightweight version 
of codegen for Python, and ``lambdarepr``, which supports many libraries 
(like NumPy), and theano (``theano_code``).The code printers 
are special cases of the other prints in SymPy (str printer, pretty printer, etc.).

An important distinction is that the code printer has to deal with 
assignments (using the ``sympy.printing.codeprinter.Assignment`` object).This serves
as building blocks for the code printers and hence the ``codegen`` module.
An example that shows the use of ``Assignment``::

    >>> from sympy import symbols, MatrixSymbol, Matrix
    >>> from sympy.printing.codeprinter import Assignment
    >>> x, y, z = symbols('x y z')
    >>> mat = Matrix([x, y, z]).T
    >>> known_mat = MatrixSymbol('K', 1, 3)
    >>> Assignment(known_mat, mat)  # doctest: +SKIP
    Assignment(K, Matrix([[x, y, z]]))
    >>> Assignment(known_mat, mat).lhs
    K
    >>> Assignment(known_mat, mat).rhs  # doctest: +SKIP
    Matrix([[x, y, z]])

Examples::

    >>> from sympy import ccode, symbols, Rational
    >>> Z, k, e, r = symbols('Z k e r')
    >>> expr = (Rational(-1, 2)*Z*k*(e**2)/r)
    >>> expr # doctest: +SKIP
        2   
    -Z⋅e ⋅k 
    ────────
      2⋅r   
    >>> ccode(expr) # doctest: +SKIP
    -1.0L/2.0L*Z*pow(e, 2)*k/r
    >>> ccode(expr, assign_to="E") # doctest: +SKIP
    E = -1.0L/2.0L*Z*pow(e, 2)*k/r;

``Piecewise`` expressions are converted into conditionals. If an
``assign_to`` variable is provided an if statement is created, otherwise
the ternary operator is used. Note that if the ``Piecewise`` lacks a
default term, represented by ``(expr, True)`` then an error will be thrown.
This is to prevent generating an expression that may not evaluate to
anything. A use case for ``Piecewise``::

    >>> from sympy import symbols, fcode, Piecewise
    >>> x, tau = symbols('x, tau')
    >>> expr = Piecewise((x + 1, x > 0), (x, True))
    >>> print(fcode(expr, tau))
          if (x > 0) then
             tau = x + 1
          else
             tau = x
          end if

The various printers also tend to support ``Indexed`` objects well.

With ``contract=True`` these expressions will be turned into loops, whereas
``contract=False`` will just print the assignment expression that should be
looped over::

    >>> from sympy import Eq, IndexedBase, Idx
    >>> from sympy import jscode
    >>> len_y = 5
    >>> y = IndexedBase('y', shape=(len_y,))
    >>> t = IndexedBase('t', shape=(len_y,))
    >>> Dy = IndexedBase('Dy', shape=(len_y-1,))
    >>> i = Idx('i', len_y-1)
    >>> e = Eq(Dy[i], (y[i+1]-y[i])/(t[i+1]-t[i]))
    >>> jscode(e.rhs, assign_to=e.lhs, contract=False) # doctest: +SKIP
    Dy[i] = (y[i + 1] - y[i])/(t[i + 1] - t[i]);

    >>> Res = IndexedBase('Res', shape=(len_y,))
    >>> j = Idx('j', len_y)
    >>> e = Eq(Res[j], y[j]*t[j])
    >>> print(jscode(e.rhs, assign_to=e.lhs, contract=True))
    for (var j=0; j<5; j++){
       Res[j] = 0;
    }
    for (var j=0; j<5; j++){
       for (var j=0; j<5; j++){
          Res[j] = Res[j] + t[j]*y[j];
       }
    }
    >>> print(jscode(e.rhs, assign_to=e.lhs, contract=False)) # doctest: +SKIP
    Res[j] = t[j]*y[j];

Custom printing can be defined for certain types by passing a dictionary of
"type" : "function" to the ``user_functions`` kwarg.  Alternatively, the
dictionary value can be a list of tuples i.e., [(argument_test,
cfunction_string)].  This can be used to call a custom Octave function::

    >>> from sympy import Function, octave_code, Function, Matrix, symbols
    >>> f = Function('f')
    >>> g = Function('g')
    >>> x = symbols('x')
    >>> custom_functions = {
    ...   "f": "existing_octave_fcn",
    ...   "g": [(lambda x: x.is_Matrix, "my_mat_fcn"),
    ...         (lambda x: not x.is_Matrix, "my_fcn")]
    ... }
    >>> mat = Matrix([[1, x]])
    >>> octave_code(f(x) + g(x) + g(mat), user_functions=custom_functions)  # doctest: +SKIP
    'existing_octave_fcn(x) + my_fcn(x) + my_mat_fcn([1 x])'


An example of mathematica code printer::

    >>> from sympy import mathematica_code as mc
    >>> from sympy import summation, symbols
    >>> from sympy import sin, Function, pprint, summation
    >>> x = Function('x')
    >>> n, T, t = symbols('n T t')
    >>> e = x(n*T)*sin((t-n*T)/T)
    >>> e = e/((-T*n + t)/T)
    >>> e   # doctest: +SKIP
    T*x(T*n)*sin((-T*n + t)/T)/(-T*n + t)
    >>> pprint(e)   # doctest: +SKIP
                ⎛-T⋅n + t⎞
    T⋅x(T⋅n)⋅sin⎜────────⎟
                ⎝   T    ⎠
    ──────────────────────
           -T⋅n + t     

    >>> expr = summation(e, (n, -1, 1))
    >>> pprint(mc(expr))    # doctest: +SKIP
    T*x[-T]*Sin[(T + t)/T]/(T + t) + T*x[T]*Sin[(-T + t)/T]/(-T + t) + T*x[0]*Sin[
    t/T]/t



We can go through a common expression in different languages we 
support and see how it works::

    >>> from sympy import jscode, ccode, fcode, octave_code, mathematica_code as mc
    >>> from sympy import cos, symbols
    >>> from sympy import pprint
    >>> k_i, gamma_i, gamma_s, r_is, I_z, S_z = symbols("k_i, gamma_i, gamma_s, r_is, I_z, S_z")
    >>> beta = symbols("beta")
    >>> e = k_i*gamma_i*gamma_s/(r_is**3)
    >>> expr = e*2*I_z*S_z*(3*(cos(beta))**2 - 1)/2
    >>> from sympy import init_printing
    >>> init_printing()
    >>> pprint(expr)    # doctest: +SKIP
                     ⎛     2       ⎞
    I_z⋅S_z⋅γᵢ⋅γₛ⋅kᵢ⋅⎝3⋅cos (β) - 1⎠
    ────────────────────────────────
                     3              
                  rᵢₛ               
    >>> pprint(jscode(expr, assign_to="H_is"))  # doctest: +SKIP
    H_is = I_z*S_z*gamma_i*gamma_s*k_i*(3*Math.pow(Math.cos(beta), 2) - 1)/Math.po
    w(r_is, 3);
    >>> pprint(ccode(expr, assign_to="H_is"))   # doctest: +SKIP
    H_is = I_z*S_z*gamma_i*gamma_s*k_i*(3*pow(cos(beta), 2) - 1)/pow(r_is, 3);
    >>> pprint(fcode(expr, assign_to="H_is"))   # doctest: +SKIP
          H_is = I_z*S_z*gamma_i*gamma_s*k_i*(3*cos(beta)**2 - 1)/r_is**3
    >>> pprint(octave_code(expr, assign_to="H_is")) # doctest: +SKIP
    H_is = I_z.*S_z.*gamma_i.*gamma_s.*k_i.*(3*cos(beta).^2 - 1)./r_is.^3;
    >>> pprint(mc(expr))    # doctest: +SKIP
    I_z*S_z*gamma_i*gamma_s*k_i*(3*Cos[beta]^2 - 1)/r_is^3

Codegen (sympy.utilities.codegen)
---------------------------------
This module deals with creating compilable code from SymPy expressions. 
This is lower level than autowrap, as it doesn't actually attempt to 
compile the code, but higher level than the printers, as it generates 
compilable files (including header files), rather than just code snippets.

The user friendly functions, here, are ``codegen`` and ``make_routine``.
``codegen`` takes a list of ``(variable, expression)`` pairs and a language 
(C, F95, and Octave/Matlab are supported). It returns, as strings, a code 
file and a header file (for relevant languages). The variables are created 
as functions that return the value of the expression as output.

.. note:: The ``codegen`` callable is not in the sympy namespace automatically,
   to use it you must first import ``codegen`` from ``sympy.utilities.codegen``

For instance::

    >>> from sympy.utilities.codegen import codegen
    >>> from sympy import symbols
    >>> length, breadth, height = symbols('length, breadth, height')
    >>> [(c_name, c_code), (h_name, c_header)] = codegen(('volume', length*breadth*height), "C", "test", header=False, empty=False)
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

Various flags to ``codegen`` let you modify things. The project name for preprocessor 
instructions can be varied using ``project``. Variables listed as global variables in 
arg ``global_vars`` will not show up as function arguments.

``language`` is a case-insensitive string that indicates the source code language. 
Currently, 'C', 'F95' and 'Octave' are supported. 
'Octave' generates code compatible with both Octave and Matlab.

``header`` when True, a header is written on top of each source file. ``empty`` 
when True, empty lines are used to structure the code. With ``argument_sequence``
a sequence of arguments for the routine can be defined in a preferred order.  

``prefix`` defines a prefix for the names of the files that contain the source code. 
If omitted, the name of the first name_expr tuple is used.
``to_files`` when True, the code will be written to one or more files with the
given prefix.
          


Here is an example::

    >>> [(f_name, f_code), header] = codegen(("volume", length*breadth*height), "F95", header=True, empty=False, argument_sequence=(breadth, length), global_vars=(height,))
    >>> print(f_code)
    !******************************************************************************
    !*                    Code generated with sympy 0.7.7.dev                     *
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



The method ``make_routine`` creates a ``Routine`` object, which represents an evaluation
routine for a set of expressions. This is only good for internal use by the CodeGen 
objects, as an intermediate representation from SymPy expression to generated code. 
It is not recommended to make a ``Routine`` object yourself. You should instead use 
``make_routine`` method. ``make_routine`` in turn calls the ``routine`` method of 
the CodeGen object depending upon the language of choice. This creates the internal 
objects representing assignments and so on, and creates the ``Routine`` class with them.

The various codegen objects such as ``Routine`` and ``Variable`` aren't SymPy 
objects (they don't subclass from Basic).

For example::

    >>> from sympy.utilities.codegen import make_routine
    >>> from sympy.physics.hydrogen import R_nl
    >>> from sympy import symbols, init_printing
    >>> init_printing()
    >>> x, y = symbols('x y')
    >>> expr = R_nl(3, y, x, 6)
    >>> r = make_routine('my_routine', expr)
    >>> [arg.result_var for arg in r.results]   # doctest: +SKIP
    [result₅₁₄₂₃₄₁₆₈₁₃₉₇₇₁₉₄₂₈]
    >>> [arg.expr for arg in r.results]
    ⎡                ___________                                           ⎤
    ⎢          y    ╱ (-y + 2)!   -2⋅x                                     ⎥
    ⎢4⋅√6⋅(4⋅x) ⋅  ╱  ───────── ⋅ℯ    ⋅assoc_laguerre(-y + 2, 2⋅y + 1, 4⋅x)⎥
    ⎢            ╲╱    (y + 3)!                                            ⎥
    ⎢──────────────────────────────────────────────────────────────────────⎥
    ⎣                                  3                                   ⎦
    >>> [arg.name for arg in r.arguments]   # doctest: +SKIP
    [x, y]

Another more complicated example with a mixture of specified and
automatically-assigned names.  Also has Matrix output::

    >>> from sympy import Matrix
    >>> from sympy.abc import x, y, f, g
    >>> r = make_routine('fcn', [x*y, Eq(f, 1), Eq(g, x + g), Matrix([[x, 2]])])
    >>> [arg.result_var for arg in r.results]   # doctest: +SKIP
    [result_5397460570204848505]
    >>> [arg.expr for arg in r.results] # doctest: +SKIP
    [x*y]
    >>> [arg.name for arg in r.arguments]   # doctest: +SKIP
    [x, y, f, g, out_8598435338387848786]

We can examine the various arguments more closely::

    >>> from sympy.utilities.codegen import (InputArgument, OutputArgument,
    ...                                      InOutArgument)
    >>> [a.name for a in r.arguments if isinstance(a, InputArgument)]   
    [x, y]

    >>> [a.name for a in r.arguments if isinstance(a, OutputArgument)]  # doctest: +SKIP
    [f, out_8598435338387848786]
    >>> [a.expr for a in r.arguments if isinstance(a, OutputArgument)]  # doctest: +SKIP
    [1, Matrix([[x, 2]])]

    >>> [a.name for a in r.arguments if isinstance(a, InOutArgument)]
    [g]
    >>> [a.expr for a in r.arguments if isinstance(a, InOutArgument)]
    [g + x]



Autowrap
--------
Autowrap automatically generates code, writes it to disk, compiles it, 
and imports it into the current session. Main functions of this module are 
``autowrap``, ``binary_function``, and ``ufuncify``.

It also automatically converts expressions containing ``Indexed`` objects 
into summations. The classes IndexedBase, Indexed and Idx represent a matrix 
element M[i, j]. See :ref:`tensor_module` for more on this.
``autowrap`` creates a wrapper using f2py or Cython and creates a numerical 
function.

.. note:: The ``autowrap`` callable is not in the sympy namespace automatically,
   to use it you must first import ``autowrap`` from ``sympy.utilities.autowrap``


The callable returned from autowrap() is a binary python function, not a 
SymPy object. For example::

    >>> from sympy.abc import x, y, z
    >>> from sympy.utilities.autowrap import autowrap
    >>> expr = ((x - y + z)**(13)).expand()
    >>> binary_func = autowrap(expr)    # doctest: +SKIP
    >>> binary_func(1, 4, 2)    # doctest: +SKIP
    -1.0

The various flags available with autowrap() help to modify the services 
provided by the method. 
The argument ‘tempdir’ tells autowrap to compile the code in a specific 
directory, and leave the files intact when finished. For instance::

    >>> from sympy.utilities.autowrap import autowrap
    >>> from sympy.physics.qho_1d import psi_n
    >>> from sympy import IndexedBase, Idx
    >>> from sympy import Eq
    >>> from sympy import symbols
    >>> x = IndexedBase('x')
    >>> y = IndexedBase('y')
    >>> m = symbols('m', integer=True)
    >>> i = Idx('i', m)
    >>> a,omega = symbols('a, omega')
    >>> qho = autowrap(Eq(y[i], psi_n(0, x[i], m, omega)), tempdir='/tmp')  # doctest: +SKIP

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

    >>> qho = autowrap(Eq(y[i], psi_n(0, x[i], m, omega)), tempdir='/tmp', args=[y, x, m, omega])   # doctest: +SKIP

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

The argument ``verbose`` is boolean, optional and if True, autowrap 
will not mute the command line backends. This can be helpful for debugging.

The argument ``language`` and ``backend`` are used to change defaults: 'Fortran'
and 'f2py' to 'C' and 'Cython'.
The argument helpers is used to define auxillary expressions needed for the main 
expression. If the main expression needs to call a specialized function it should 
be put in the ``helpers`` iterable. Autowrap will then make sure that the
compiled main expression can link to the helper routine. Items should
be tuples with (<function_name>, <sympy_expression>, <arguments>). It is mandatory 
to supply an argument sequence to helper routines.

Another method available at the ``autowrap`` level is ``binary_function``. It returns 
a sympy function. The advantage is that we can have very fast functions as compared
to SymPy speeds. This is because we will be using compiled functions with Sympy attriutes 
and methods. An illustration::

    >>> from sympy.utilities.autowrap import binary_function
    >>> from sympy import symbols
    >>> from sympy.physics.hydrogen import R_nl
    >>> a, r = symbols('a, r')
    >>> psi_nl = R_nl(1, 0, a, r)
    >>> f = binary_function('f', psi_nl)    # doctest: +SKIP
    >>> f(a, r).evalf(3, subs={a: 1, r: 2})  # doctest: +SKIP
    0.766


While NumPy operations are very efficient for vectorized data but they sometimes incur 
unnecessary costs when `chained together <http://docs.sympy.org/dev/modules/numeric-computation.html#ufuncify>`_  .
Fortunately, SymPy is able to generate efficient low-level C or Fortran code. 
It can then depend on projects like Cython or f2py to compile and reconnect that 
code back up to Python. Fortunately this process is well automated and a SymPy user 
wishing to make use of this code generation should call the ufuncify function.
``ufuncify`` is the third method available with Autowrap module. 
It basically implies 'Universal functions' and follows an ideology set by Numpy.
The main point of ufuncify as compared to autowrap is that it allows arrays as arguments 
and can operate in an element-by-element fashion. The core operation done element-wise is 
in accordance to Numpy's array broadcasting rules.
See `this <http://docs.scipy.org/doc/numpy/reference/ufuncs.html>`_ for more.

Let us see an example::

    >>> from sympy import init_printing, symbols
    >>> init_printing()
    >>> from sympy.physics.hydrogen import R_nl
    >>> x = symbols('x')
    >>> expr = R_nl(3, 1, x, 6)
    >>> expr
                    -2⋅x
    8⋅x⋅(-4⋅x + 4)⋅ℯ    
    ────────────────────
             3          


The lambdify function translates SymPy expressions into Python functions, 
leveraging a variety of numerical libraries.By default lambdify relies 
on implementations in the ``math`` standard library. Naturally, Raw Python 
is faster than Sympy. However it also supports ``mpmath`` and most notably, 
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