=================
SymPy vs. Yacas
=================

|SymPy| vs. |Yacas|

.. |SymPy| image:: SymPy.png
.. |Yacas| image:: Yacas.png

SymPy_ and Yacas_ are *Computer algebra systems*.

.. _SymPy: http://sympy.org/
.. _Yacas: http://yacas.sourceforge.net/homepage.html

**Computer Algebra System**
    A software program that facilitates symbolic mathematics.
    The core functionality of a CAS is manipulation of mathematical expressions in symbolic form.

+++++++
SymPy
+++++++

**Sympy** is a Python library for symbolic computation that aims to become a full-featured computer algebra system and to keep the code simple to promote extensibility and comprehensibility.

SymPy was started by Ondřej Čertík in 2005 and he wrote some code in 2006 as well. In 11 March 2007, SymPy was realeased to the public.
The latest stable release of SymPy is 0.7.1 (29 July 2011). As of beginning of December 2011 there have been over 150 people who contributed at least one commit to SymPy.

SymPy can be used:

- Inside Python, as a library
- As an interactive command line, using IPython

SymPy is entirely written in Python and does not require any external libraries, but various programs that can extend its capabilites can be installed:

- gmpy, Cython --> speed improvement
- Pyglet, Matplotlib --> 2d and 3d plotting
- IPython --> interactive sessions

SymPy is available online at `SymPy Live`_. The site was developed specifically for SymPy. It is a simple web shell that looks similar to iSymPy under the standard Python interpreter. SymPy Live uses Google App Engine as computational backend.

.. _`SymPy Live`: http://live.sympy.org/

\+ \+: small library, pure Python, very functional, extensible, large community.

\- \-: slow, needs better documentation.

+++++++
Yacas
+++++++

**Yacas** is a small and highly flexible general-purpose computer algebra system and programming language. The name is an acronym for Yet Another Computer Algebra System. Yacas is a program for symbolic manipulation of mathematical expressions.

It uses its own programming language, designed for symbolic as well as arbitrary-precision numerical computations. The system has a library of scripts that implement many of the symbolic algebra operations. New algorithms can be easily added to the library.

Yacas was created by Ayal Pinkus et al. in 1998. The latest stable release of Yacas is 1.3.1 (25 October 2011). Yacas is written entirely in C++.

Yacas is primarily intended to be a research tool for easy exploration and prototyping of algorithms of symbolic computation. The main advantage of Yacas is its rich and flexible scripting language.

It is possible to try Yacas online at the `My Yacas`_ section of Yacas' official site.

.. _`My Yacas`: http://yacas.sourceforge.net/homepage.html

\+ \+: small library, very functional, flexible.

\- \-: slow, needs better documentation, small community.

++++++++++++++++++
Sympy (!)= Yacas
++++++++++++++++++

Both *SymPy* and *Yacas* are cost free open source CASes. SymPy is released under a modified BSD license, while *Yacas* is realeased unde the GNU GPL.

One of the differences between SymPy and Yacas is that Yacas has a GUI for Windows platforms called `GUYacas`_. This program is selfcontainted, meaning that it isn't necessary to install Yacas for GUYacas to work. However, GUYacas requires that 'Microsoft .NET Framework Version 2.0' is installed.

.. _`GUYacas`: http://gbi.agrsci.dk/~sorenh/misc/GUYacas/index.html

-------------------------
Operating System Support
-------------------------

+------------+---------+----------+-------+-----+---------+-----------------------------------+
| System     | Windows | Mac OS X | Linux | BSD | Solaris |               Other               |
+------------+---------+----------+-------+-----+---------+-----------------------------------+
|  SymPy     |   Yes   |    Yes   |  Yes  | Yes |   Yes   |  Any system that supports Python  |
+------------+---------+----------+-------+-----+---------+-----------------------------------+
|  Yacas     |   Yes   |    Yes   |  Yes  | Yes |   Yes   |                 ?                 |
+------------+---------+----------+-------+-----+---------+-----------------------------------+

------------------------
Download & Installation
------------------------

Sympy is distributed in various forms. It is possible to download source tarballs and packages from the Google Code page but it is also possible to clone the main Git repository or browse the code online. The only prerequisite is Python since Sympy is Python-based library. It is recommended to install IPython as well, for a better experience.

Yacas is small (1.2 MB) and portable across most Unix-like platforms and requires a standard C++ compiler such as g++. Precompiled Red Hat (RPM) and Debian (DEB) packages are also available. To compile and install Yacas, read the file INSTALL. There are two methods to compile and install Yacas: downloading a tarball or downloading from CVS at SourceForge.

--------------
Functionality
--------------

+------------+----------+------------+-----------------------------------+---------------------------------------------------------------------------+
|            | Formula  | Arbitrary  |             Calculus              |                                            Solvers                        |
|  System    |          |            +-------------+---------------------+-----------+--------------+-----------------------+------------------------+
|            | editor   | precision  | Integration |Integral transforms* | Equations | Inequalities | Diophantine equations | Differential equations |
+------------+----------+------------+-------------+---------------------+-----------+--------------+-----------------------+------------------------+
|  SymPy     |    No    |    Yes     |    Yes      |        No           |   Yes     |     Yes      |          No           |           Yes          |
+------------+----------+------------+-------------+---------------------+-----------+--------------+-----------------------+------------------------+
|  Yacas     |    No    |    Yes     |    Yes      |        No           |   Yes     |     No       |          No           |           No           |
+------------+----------+------------+-------------+---------------------+-----------+--------------+-----------------------+------------------------+

+------------+-----------------------+---------+---------+--------------+----------+---------+
|            |        Solvers        | Graph   | Number  | Quantifier   | Boolean  |         |
|  System    +-----------------------+         |         |              |          | Tensors |
|            | Recurrence relations  | theory  | theory  | elimination  | algebra  |         |
+------------+-----------------------+---------+---------+--------------+----------+---------+
|  SymPy     |          Yes          |   No    |   Yes   |     No       |   Yes    |   Yes   |
+------------+-----------------------+---------+---------+--------------+----------+---------+
|  Yacas     |          No           |   No    |   No    |     No       |   No     |   No    |
+------------+-----------------------+---------+---------+--------------+----------+---------+

\* Will be available in SymPy 0.7.2

''''''''''''''''''''''''''
Some syntax differences
''''''''''''''''''''''''''

In SymPy, to raise something to a power, you must use \*\*, not ^ as the latter uses the Python meaning, which is xor.

::

    In [1]: (x+1)^2
    ---------------------------------------------------------------------------
    TypeError                                 Traceback (most recent call last)
    /home/aoi_hana/sympy/<ipython-input-6-52730bce1577> in <module>()
    ----> 1 (x+1)^2

    TypeError: unsupported operand type(s) for ^: 'Add' and 'int'

    In [2]: (x+1)**2
    Out[2]:
           2
    (x + 1)

However, in Yacas, you must use ^ for exponentiation, as \*\* isn't defined as an operator:

::

    In> (x+1)**2
    String(1) : Error parsing expression, near token **

    In> (x+1)^2
    Out> (x+1)^2

Another difference between SymPy and Yacas is that you have to define symbols in SymPy before you can use them, while in Yacas it isn't necessary.

**SymPy**

::

    >>> x**2 + 2*x + 1
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    NameError: name 'x' is not defined

    >>> from sympy import Symbol
    >>> x = Symbol('x')
    >>> x**2 + 2*x + 1
    x**2 + 2*x + 1

**Yacas**

    In> x^2+2*x+1
    Out> x^2+2*x+1

''''''''''
Algebra
''''''''''

**SymPy**

To perform partial fraction decomposition *apart(expr, x)* must be used. To combine expressions, *together(expr, x)* is what you need.
Here are some examples of these two and other common functions in iSymPy:

::

    In [8]: 1/( (x**2+2*x+1)*(x**2-1) )
    Out[8]:
               1
    ───────────────────────
    ⎛ 2    ⎞ ⎛ 2          ⎞
    ⎝x  - 1⎠⋅⎝x  + 2⋅x + 1⎠

    In [9]: apart(1/( (x**2+2*x+1)*(x**2-1) ), x)
    Out[9]:
          1           1            1            1
    - ───────── - ────────── - ────────── + ─────────
      8⋅(x + 1)            2            3   8⋅(x - 1)
                  4⋅(x + 1)    2⋅(x + 1)

    In [10]: together(1/(x**2+2*x) - 3/(x+y) + 1/(x+y+z))
    Out[10]:
    x⋅(x + 2)⋅(x + y) - 3⋅x⋅(x + 2)⋅(x + y + z) + (x + y)⋅(x + y + z)
    ─────────────────────────────────────────────────────────────────
                      x⋅(x + 2)⋅(x + y)⋅(x + y + z)

The *evalf()* method and the *N()* function can be used to evaluate expressions:

::

    In [20]: pi.evalf()
    Out[20]: 3.14159265358979

    In [23]: N(sqrt(2)*pi, 50)
    Out[23]: 4.4428829381583662470158809900606936986146216893757

Integrals can be used like regular expressions and support arbitrary precision:

::

    In [24]: Integral(x**(-2*x), (x, 0, oo)).evalf(20)
    Out[24]: 2.0784499818221828310

**Yacas**

The *Apart(expr)* and *Together(expr)* functions work in Yacas as well.

::

    In> PrettyPrinter'Set("PrettyForm")

    True

    In> 1/( (x^2+2*x+1)*(x^2-1) )

                   1
    -------------------------------
    /  2             \   /  2     \
    \ x  + 2 * x + 1 / * \ x  - 1 /

    In> Apart(1/( (x^2+2*x+1)*(x^2-1) ))

          1         /       1                1                1       \
    ------------- - | -------------- + -------------- + ------------- |
    8 * ( x - 1 )   |              3                2   8 * ( x + 1 ) |
                    \ 2 * ( x + 1 )    4 * ( x + 1 )                  /

    In> Together(1/(x^2+2*x)- 3/(x+y) + 1/(x+y+z))

    /     1          3   \
    | ---------- - ----- | * ( x + y + z ) + 1
    |  2           x + y |
    \ x  + 2 * x         /
    ------------------------------------------
                    x + y + z

*N(expr)* or *N(expr, prec)* can be used to determine the numerical approximation of expressions.

::

    In> N(Pi)
    Out> 3.1415926535897932384

    In> N(Sqrt(2)*Pi, 49)
    Out> 4.4428829381583662470158809900606936986146216893757

''''''''''
Calculus
''''''''''

""""""""""
Limits
""""""""""

**SymPy**

Limits in SymPy have the following syntax: *limit(function, variable, point)*.
Here are some examples:

Limit of f(x)= sin(x)/x as x -> 0

::

    In [20]: from sympy import *

    In [21]: x = Symbol('x')

    In [22]: limit(sin(x)/x, x, 0)
    Out[22]: 1

Limit of f(x)= 2*x+1 as x -> 5/2

::

    In [24]: limit(2*x+1, x, S(5)/2)     # The *S()* method must be used for 5/2 to be Rational in SymPy
    Out[24]: 6

You can also compute the left and right limits of an expression with the *dir="+/-"* argument.

::

    In [5]: limit(1/x, x, oo)
    Out[5]: 0

    In [6]: limit(1/x, x, 0, dir="+")
    Out[6]: ∞

    In [7]: limit(1/x, x, 0, dir="-")
    Out[7]: -∞

**Yacas**

The function *Limit(var, val) expr* returns the limit of an expression.

::

    In> Limit(x,0) Sin(x)/x
    Out> 1

    In> Limit(x,5/2) 2*x+1
    Out> 6

It is also possible to compute the left and right limits of a given expression.

::

    In> Limit(x,0) 1/x
    Out> Undefined
    In> Limit(x,0,Left) 1/x
    Out> -Infinity
    In> Limit(x,0,Right) 1/x
    Out> Infinity

"""""""""""""""""
Differentiation
"""""""""""""""""

**SymPy**

::

    In [1]: from sympy import *

    In [2]: x = Symbol('x')

    In [3]: diff(cos(x**3), x)
    Out[3]:
        2    ⎛ 3⎞
    -3⋅x ⋅sin⎝x ⎠

    In [4]: diff(atan(2*x), x)
    Out[4]:
       2
    ────────
       2
    4⋅x  + 1

    In [6]: diff(1/tan(x), x)
    Out[6]:
         2
    - tan (x) - 1
    ─────────────
         2
      tan (x)

This is how you create a Bessel function of the first kind object and differentiate it:

::

    In [7]: from sympy import besselj, jn

    In [8]: from sympy.abc import z, n

    In [9]: b = besselj(n, z)

    In [10]: # Differentiate it:

    In [11]: b.diff(z)
    Out[11]:
    besselj(n - 1, z)   besselj(n + 1, z)
    ───────────────── - ─────────────────
            2                   2

**Yacas**

*D(var) expr* calculates the derivative of the expression 'expr' with respect to the variable 'var' and returns it.

::

    In> D(x) Cos(x^3)
    Out> -3*x^2*Sin(x^3)

    In> PrettyForm(D(x) ArcTan(2*x))

          2
    --------------
             2
    ( 2 * x )  + 1

    Out> True

    In> PrettyForm(D(x) 1/Tan(x))

             -1
    ---------------------
            2           2
    Cos( x )  * Tan( x )

    Out> True

""""""""""""""""""
Series expansion
""""""""""""""""""

**SymPy**

The syntax for series expansion is: *.series(var, point, order)*:

::

    In [27]: from sympy import *

    In [28]: x = Symbol('x')

    In [29]: cos(x).series(x, 0, 14)
    Out[29]:
         2    4     6      8       10         12
        x    x     x      x       x          x         ⎛ 14⎞
    1 - ── + ── - ─── + ───── - ─────── + ───────── + O⎝x  ⎠
        2    24   720   40320   3628800   479001600

    In [30]: (1/cos(x**2)).series(x, 0, 14)
    Out[30]:
         4      8       12
        x    5⋅x    61⋅x      ⎛ 14⎞
    1 + ── + ──── + ────── + O⎝x  ⎠
        2     24     720

It is possible to make use of *series(x*cos(x), x)* by creating a wrapper around Basic.series().

::

    In [31]: from sympy import Symbol, cos, series
    In [32]: x = Symbol('x')
    In [33]: series(cos(x), x)
    Out[33]:
         2    4
        x    x     ⎛ 6⎞
    1 - ── + ── + O⎝x ⎠
        2    24

This module also implements automatic keeping track of the order of your expansion.

::

    In [1]: from sympy import Symbol, Order

    In [2]: x = Symbol('x')

    In [3]: Order(x) + x**2
    Out[3]: O(x)

    In [4]: Order(x) + 28
    Out[4]: 28 + O(x)

**Yacas**

*Taylor(var, at, order) expr* returns the Taylor series expansion of the expression 'expr' with respect to the variable 'var' around 'at' up to order 'order'.

::

    In> PrettyForm(Taylor(x,0,14) Cos(x))

         2    4    6      8        10         12          14
        x    x    x      x        x          x           x
    1 - -- + -- - --- + ----- - ------- + --------- - -----------
        2    24   720   40320   3628800   479001600   87178291200

    Out> True

    In> PrettyForm(Taylor(x,0,14) Sin(x))

         3    5      7       9        11          13
        x    x      x       x        x           x
    x - -- + --- - ---- + ------ - -------- + ----------
        6    120   5040   362880   39916800   6227020800

"""""""""""""
Integration
"""""""""""""

**SymPy**

The *integrals* module in SymPy implements methods calculating definite and indefinite integrals of expressions.
Principal method in this module is *integrate()*:

- integrate(f, x) returns the indefinite integral |int1|
- integrate(f, (x, a, b)) returns the definite integral |int2|

.. |int1| image:: int1.png
.. |int2| image:: int2.png

SymPy can integrate:

- polynomial functions:

::

    In [6]: from sympy import *

    In [7]: import sys

    In [8]: from sympy import init_printing

    In [9]: init_printing(use_unicode=False, wrap_line=False, no_global=True)

    In [10]: x = Symbol('x')

    In [11]: integrate(x**2 + 2*x + 4, x)
     3
    x     2
    ── + x  + 4⋅x
    3

- rational functions:

::

    In [1]: integrate((x+1)/(x**2+4*x+4), x)
    Out[1]:
                   1
    log(x + 2) + ─────
                 x + 2

- exponential-polynomial functions:

::

    In [5]: integrate(5*x**2 * exp(x) * sin(x), x)
    Out[5]:
       2  x             2  x                             x             x
    5⋅x ⋅ℯ ⋅sin(x)   5⋅x ⋅ℯ ⋅cos(x)        x          5⋅ℯ ⋅sin(x)   5⋅ℯ ⋅cos(x)
    ────────────── - ────────────── + 5⋅x⋅ℯ ⋅cos(x) - ─────────── - ──────────
          2                2                               2             2

- non-elementary integrals:

::

    In [11]: integrate(exp(-x**2)*erf(x), x)
      ___    2
    ╲╱ π ⋅erf (x)
    ─────────────
          4

Here is an example of a definite integral (Calculate |integral1|):

.. |integral1| image:: int3.png

::

    In [1]: integrate(x**2 * cos(x), (x, 0, pi/2))
    Out[1]:
          2
         π
    -2 + ──
         4

**Yacas**

*Integrate(var, x1, x2) expr* and *Integrate(var) expr* are the for definite and indefinite integration. Here are a few examples:

- polynomial functions:

::

    In> Integrate(x) x^2+2*x+4
    Out> x^3/3+x^2+4*x

- rational functions:

::

    In> PrettyForm(Integrate(x) (x+1)/(x^2+4*x+4))

                           -1
    Ln( x + 2 ) + ( x + 2 )

    Out> True

- exponential-polynomial functions:

::

    In> Integrate(x) 5*x^2*Exp(x)*Sin(x)
    Out> AntiDeriv(x,5*Exp(x)*x^2*Sin(x))

    In> Integrate(x) Exp(x)+ x^2*Sin(x)
    Out> Exp(x)-(x^2*Cos(x)+(-2)*Cos(x)-2*x*Sin(x))

- non-elementary integrals:

::

    In> 2*Sqrt(Pi)*Integrate(x,0,1) Exp(-t^2)
    Out> 2*Sqrt(Pi)*Exp(-t^2)

An examples of a definite integral:

::

    In> PrettyForm(Integrate(x,0,Pi/2) x^2*Cos(x))

          2
    / Pi \
    | -- |
    \ 2  /  - 2

    Out> True

"""""""""""""""""
Complex numbers
"""""""""""""""""

**SymPy**

::

    In [1]: from sympy import Symbol, exp, I

    In [2]: x = Symbol("x")

    In [3]: exp(I*2*x).expand()
    Out[3]:
     2⋅ⅈ⋅x
    ℯ

    In [4]: exp(I*2*x).expand(complex=True)
    Out[4]:
       -2⋅im(x)                 -2⋅im(x)
    ⅈ⋅ℯ        ⋅sin(2⋅re(x)) + ℯ        ⋅cos(2⋅re(x))

    In [5]: x = Symbol("x", real=True)

    In [6]: exp(I*2*x).expand(complex=True)
    Out[6]: ⅈ⋅sin(2⋅x) + cos(2⋅x)

    In [7]: exp(-2 + 3*I*x).expand(complex=True)
    Out[7]:
      -2             -2
    ⅈ⋅ℯ  ⋅sin(3⋅x) + ℯ  ⋅cos(3⋅x)

Complex number division in iSymPy:

::

    In [4]: from sympy import I
    In [5]: ((2 + 3*I)/(3 + 7*I)).expand(complex=True)
    Out[5]:
    27   5⋅ⅈ
    ── - ───
    58    58

**Yacas**

::

    In> Exp(I*2*x)
    Out> Exp(Complex(0,2)*x)

Complex number division in Yacas:

::

    In> PrettyForm(Complex(2,3)/Complex(3,7))

           / 27   -5 \
    Complex| -- , -- |
           \ 58   58 /

    Out> True

"""""""""""
Functions
"""""""""""

**SymPy**

**trigonometric**

::

    In [1]: cos(x-y).expand(trig=True)
    Out[1]: sin(x)⋅sin(y) + cos(x)⋅cos(y)

    In [2]: cos(2*x).expand(trig=True)
    Out[2]:
         2
    2⋅cos (x) - 1

    In [3]: sinh(I*x**2)
    Out[3]:
         ⎛ 2⎞
    ⅈ⋅sin⎝x ⎠

    In [11]: sinh(acosh(x))
    Out[11]:
      _______   _______
    ╲╱ x - 1 ⋅╲╱ x + 1

**zeta function**

::

    In [4]: zeta(5, x**2)
    Out[4]:
     ⎛    2⎞
    ζ⎝5, x ⎠

    In [5]: zeta(5, 2)
    Out[5]: ζ(5, 2)

    In [6]: zeta(4, 1)
    Out[6]:
     4
    π
    ──
    90

    In [5]: zeta(28).evalf()
    Out[5]: 1.00000000372533

**factorials and gamma function**

::

    In [7]: a = Symbol('a')

    In [8]: b = Symbol('b', integer=True)

    In [9]: factorial(a)
    Out[9]: a!

    In [10]: N(gamma(S(25)/10), 31)
    Out[10]: 1.329340388179137020473625612506

**polynomials**

::

    In [14]: chebyshevt(8,x)
    Out[14]:
         8        6        4       2
    128⋅x  - 256⋅x  + 160⋅x  - 32⋅x  + 1

    In [15]: legendre(3, x)
    Out[15]:
       3
    5⋅x    3⋅x
    ──── - ───
     2      2

    In [16]: hermite(3, x**2)
    Out[16]:
       6       2
    8⋅x  - 12⋅x

**Yacas**

**trigonometric**

*TrigSimpCombine(expr)* applies the product rules of trigonometry.

::

    In> TrigSimpCombine(Cos(x)^2-Sin(x)^2)
    Out> Cos((-2)*x)


    In> PrettyForm(TrigSimpCombine(Cos(a)^2*Sin(b)))

    Sin( b )   Sin( -2 * a + b )   Sin( -2 * a - b )
    -------- + ----------------- - -----------------
      2               4                   4

    Out> True

    In> ArcCos(Sqrt(3)/2)
    Out> Pi/6

    In> Sin(ArcSin(alpha))+Tan(ArcTan(beta))
    Out> alpha+beta

**zeta function**

*Zeta(x)* is an interface to Riemann's Zeta function zeta(s).

::

    In> PrettyForm(Zeta(5))

    Zeta( 5 )

    Out> True

    In> N(Zeta(5.2))
    Out> 1.031659876678

    In> Zeta(4)
    Out> Pi^4/90

    In> N(Zeta(28), 13)
    Out> 1.0000000037253

**factorials and gamma function**

::

    In> a!
    Out> a!

*Gamma(x)* in an interface to Euler's Gamma function Gamma(x).

::

    In> N(Gamma(2.5), 30)
    Out> 1.329340388179137020473625612505

**polynomials**

*OrthoT(n, x)* evaluates the Chebyshev polynomials of the first kind T(n,x), of degree n at the point x.

::

    In> PrettyForm(OrthoT(8, x))

    / / /        2       \    2       \    2      \    2
    \ \ \ 128 * x  - 256 / * x  + 160 / * x  - 32 / * x  + 1

    Out> True

*OrthoP(n,x)* evaluates the Legendre polynomial of degree n at the point x.

::

    In> PrettyForm(OrthoP(3,x))

        /      2     \
        | 5 * x    3 |
    x * | ------ - - |
        \   2      2 /

    Out> True

*OrthoH(n, x)* evaluates the Hermite polynomial of degree n at the point x.

::

    In> PrettyForm(OrthoH(3, x^2))

     2   /      4      \
    x  * \ 8 * x  - 12 /

    Out> True

""""""""""""""""""""""""
Differential equations
""""""""""""""""""""""""

**SymPy**

In *iSymPy*:

::

    In [10]: f(x).diff(x, x) + f(x)
    Out[10]:
             2
            d
    f(x) + ───(f(x))
             2
           dx

    In [11]: dsolve(f(x).diff(x, x) + f(x), f(x))
    Out[11]: f(x) = C₁⋅sin(x) + C₂⋅cos(x)

**Yacas**

*OdeSolve(expr1==expr2)* can solve second order homogeneous linear real constant coefficient equations.

::

    In> OdeSolve(y'+y==0)
    Out> C82*Exp(-x)
    In> OdeSolve(y''+4*y'+2*y==0)
    Out> C113*Exp((-x*(Sqrt(8)+4))/2)+C117*Exp((x*(Sqrt(8)-4))/2)

*OdeOrder(eqn)* returns the order of the differential equation, which is the order of the highest derivative. This function returns zero whne no derivatives appear.

::

    In> OdeOrder(Cos(x)*y(5)+4*y''+5*y'==0)
    Out> 5

"""""""""""""""""""""
Algebraic equations
"""""""""""""""""""""

**SymPy**

In *iSymPy*:

::

    In [3]: solve(x**3 + 2*x**2 - 1, x)
    Out[3]:
    ⎡            ___      ___    ⎤
    ⎢      1   ╲╱ 5     ╲╱ 5    1⎥
    ⎢-1, - ─ + ─────, - ───── - ─⎥
    ⎣      2     2        2     2⎦


    In [5]: solve( [x**2 + 4*y**2 -2, -10*x + 2*y -15], [x, y])
    Out[5]:
    ⎡⎛          ____              ____  ⎞  ⎛          ____              ____   ⎞⎤
    ⎢⎜  150   ╲╱ 23 ⋅ⅈ   15   5⋅╲╱ 23 ⋅ⅈ ⎟  ⎜  150   ╲╱ 23 ⋅ⅈ   15   5⋅╲╱ 23 ⋅  ⎟⎥
    ⎢⎜- ─── - ────────, ─── - ──────────⎟, ⎜- ─── + ────────, ─── + ────────── ⎟⎥
    ⎣⎝  101     101     202      101    ⎠  ⎝  101     101     202      101     ⎠⎦

**Yacas**

::

    In> Solve(2*x^2 - 1, x)

    /                     \
    |          / 1 \      |
    | x == Sqrt| - |      |
    |          \ 2 /      |
    |                     |
    |      -( Sqrt( 8 ) ) |
    | x == -------------- |
    |            4        |
    \                     /

    In> Solve({x^2+4*y^2==2, -10*x+2*y==15}, {x, y})

    /                                         \
    | /          /          2 \ \ ( y == y )  |
    | \ x == Sqrt\ 2 - 4 * y  / /             |
    \                                         /

''''''''''''''''
Linear Algebra
''''''''''''''''

""""""""""
Matrices
""""""""""

**SymPy**

In SymPy, matrices are created as instances from the Matrix class:

::

    In [1]: from sympy import Matrix

    In [2]: Matrix([ [1, 0 , 0], [0, 1, 0], [0, 0, 1] ])
    Out[2]:
    ⎡1  0  0⎤
    ⎢       ⎥
    ⎢0  1  0⎥
    ⎢       ⎥
    ⎣0  0  1⎦

It is possible to slice submatrices, since this is Python:

::

    In [4]: M = Matrix(2, 3, [1, 2, 3, 4, 5, 6])

    In [5]: M[0:2,0:2]
    Out[5]:
    ⎡1  2⎤
    ⎢    ⎥
    ⎣4  5⎦

    In [6]: M[1:2,2]
    Out[6]: [6]

    In [7]: M[:,2]
    Out[7]:
    ⎡3⎤
    ⎢ ⎥
    ⎣6⎦

One basic operation involving matrices is the determinant:

::

    In [8]: M = Matrix(( [2, 5, 6], [4, 7, 10], [1, 0, 3] ))

    In [9]: M.det()
    Out[9]: -10

*print_nonzero(symb='x')* shows location of non-zero entries for fast shape lookup.

::

    In [10]: M = Matrix(( [2, 0, 0, 1, 0], [3, 5, 0, 1, 0], [10, 4, 0, 1, 2], [1, 6, 0, 0, 0], [0, 4, 0, 2, 2] ))
    In [12]: M
    Out[12]:
    ⎡2   0  0  1  0⎤
    ⎢              ⎥
    ⎢3   5  0  1  0⎥
    ⎢              ⎥
    ⎢10  4  0  1  2⎥
    ⎢              ⎥
    ⎢1   6  0  0  0⎥
    ⎢              ⎥
    ⎣0   4  0  2  2⎦

    In [13]: M.print_nonzero()
    [X  X ]
    [XX X ]
    [XX XX]
    [XX   ]
    [ X XX]

Matrix transposition with **transpose()**:

::

    In [14]: from sympy import Matrix, I

    In [15]: m = Matrix(( (1,2+I), (3,4) ))

    In [16]: m
    Out[16]:
    ⎡1  2 + ⅈ⎤
    ⎢        ⎥
    ⎣3    4  ⎦

    In [17]: m.transpose()
    Out[17]:
    ⎡  1    3⎤
    ⎢        ⎥
    ⎣2 + ⅈ  4⎦

    In [19]: m.T == m.transpose()
    Out[19]: True

The *multiply_elementwise(b)* method returns the Hadamard product (elementwise product) of A and B:

::

    In [14]: import sympy

    In [15]: A = sympy.Matrix([ [1, 3, 20], [1, 18, 3] ])
    In [17]: B = sympy.Matrix([ [0, 5, 10], [4, 20, 6] ])

    In [18]: print A.multiply_elementwise(B)
    [0,  15, 200]
    [4, 360,  18]

**Yacas**

To make an identity matrix in Yacas, *Identity(n)* must be used, where n is the size of the matrix.

::

    In> Identity(3)
    Out> {{1,0,0},{0,1,0},{0,0,1}}

The command *DiagonalMatrix(d)* constructs a diagonal matrix, a square matrix whose off-diagonal entries are all zero.

::

    In> DiagonalMatrix(1 .. 6)
    Out> {{1,0,0,0,0,0},{0,2,0,0,0,0},{0,0,3,0,0,0},{0,0,0,4,0,0},{0,0,0,0,5,0},{0,0,0,0,0,6}}

*Determinant(M)* returns the determinant of a matrix M.

::

    In> M:={{2,5,6},{4,7,10},{1,0,3}}
    Out> {{2,5,6},{4,7,10},{1,0,3}}
    In> Determinant(M)
    Out> -10

*Transpose(M)* returns the transpose of a matrix M. This operation is useful for lists too, since matrices are just lists of lists.

::

    In> A:={{1,2+I},{3,4}}
    Out> {{1,Complex(2,1)},{3,4}}
    In> Transpose(A)
    Out> {{1,3},{Complex(2,1),4}}
    In> PrettyForm(Transpose(A))

    /                             \
    | ( 1 )                ( 3 )  |
    |                             |
    | ( Complex( 2 , 1 ) ) ( 4 )  |
    \                             /

*IsSquareMatrix(expr)* returns True if expr is a square matrix, False otherwise.

::

    In> IsSquareMatrix({{13,15},{28,82}});
    Out> True
    In> IsSquareMatrix({{13,14,15},{28,82,41}});
    Out> False

''''''''''
Geometry
''''''''''

**SymPy**

The geometry module can be used to create two-dimensional geometrical entities and query information about them.
These entities are available:

- Point
- Line, Ray, Segment
- Ellipse, Circle
- Polygon, RegularPolygon, Triangle

Check if points are collinear:

::

    In [37]: from sympy import *

    In [38]: from sympy.geometry import *

    In [39]: x = Point(0, 0)

    In [40]: y = Point(3, 1)

    In [41]: z = Point(5, 5)

    In [42]: Point.is_collinear(x, y, z)
    Out[42]: False

    In [43]: Point.is_collinear(x, z)
    Out[43]: True

Segment declaration, slope, length, midpoint:

::

    In [1]: import sympy

    In [2]: from sympy import Point

    In [3]: from sympy.abc import s

    In [4]: from sympy.geometry import Segment

    In [5]: Segment( (1, 2), (2, -3))
    Out[5]: ((1,), (2,))

    In [6]: s = Segment(Point(4, 3), Point(1, 1))

    In [7]: s
    Out[7]: ((1,), (4,))

    In [8]: s.points
    Out[8]: ((1,), (4,))

    In [9]: s.slope
    Out[9]: 2/3

    In [10]: s.length
    Out[10]:
      ____
    ╲╱ 13

    In [11]: s.midpoint
    Out[11]: (5/2,)

**Yacas**

Yacas doesn't have support for geometry yet.

''''''''''''''''''
Pattern matching
''''''''''''''''''

**SymPy**

Using the *.match* method and the *Wild* class you can perform pattern matching on expressions.
The method returns a dictionary with the needed substitutions. Here is an example:

::

    In [11]: from sympy import *

    In [12]: x = Symbol('x')

    In [13]: y = Wild('y')

    In [14]: (10*x**3).match(y*x**3)
    Out[14]: {y: 10}

    In [15]: s = Wild('s')

    In [16]: (x**4).match(y*x**s)
    Out[16]: {s: 4, y: 1}

SymPy returns *None* if the match is unsuccessful:

::

    In [19]: print (x+1).match(y**x)
    None

**Yacas**

::

    In> f(0) <-- 1

    True

    In> f(n_IsPositiveInteger) <-- n*f(n-1)

    True

    In> log(_x * _y) <-- log(x) + log(y)

    True

    In> log(_x ^ _n) <-- n * log(x)

    True

*MatchLinear(x, expr)* tries to match an expression to a linear polynomial.

::

    In> MatchLinear(x, (R+1)*x+(T-1))

    True

    In> MatchLinear(x, Sin(x)*x+(T-1))

    False

''''''''''
Printing
''''''''''

**SymPy**

There are many ways of printing mathematical expressions.
Three of the most common methods are:

- Standard printing
- Pretty printing using the pprint() function
- Pretty printing using the init_printing() method

*Standard printing* is the return value of *str(expression)*:

::

    >>> from sympy import Integral   # Python session
    >>> from sympy.abc import c
    >>> print c**3
    c**3
    >>> print 2/c
    2/c
    >>> print Integral(c**2+2*c, c)
    Integral(c**2 + 2*c, c)

*Pretty printing* is a nice ascii-art printing with the help of a *pprint* function.

::

    In [1]: from sympy import Integral, pprint   # IPython session (pprint enabled by default)

    In [2]: from sympy.abc import c

    In [3]: pprint(c**3)
     3
    c

    In [4]: pprint(2/c)
    2
    ─
    c

    In [5]: pprint(Integral(c**2+2*c, c))
    ⌠
    ⎮ ⎛ 2      ⎞
    ⎮ ⎝c  + 2⋅c⎠ dc
    ⌡

However, the proper way to set up pretty printing in SymPy is to use *init_printing(pretty_print=True, order=None, use_unicode=None, wrap_line=None, num_columns=None, no_global=False, ip=None)*:

::

    >>> from sympy import init_printing
    >>> init_printing(use_unicode=False, wrap_line=False, no_global=True)
    >>> from sympy import Integral, Symbol
    >>> x = Symbol('x')
    >>> Integral(x**3+2*x+1, x)
      /
     |
     | / 3          \
     | \x  + 2*x + 1/ dx
     |
    /
    >>> init_printing(pretty_print=True)
    >>> Integral(x**3+2*x+1, x)
    ⌠
    ⎮ ⎛ 3          ⎞
    ⎮ ⎝x  + 2⋅x + 1⎠ dx
    ⌡

**Yacas**

There are several methods to print expressions in Yacas.
The most common methods are listed below:

- FullForm - print an expression in LISP-format
- Echo - high-level printing routine
- PrettyForm - print an expression nicely with ASCII art
- EvalFormula - print an evaluation nicely with ASCII art
- TeXForm - export expressions to LaTex

*FullForm(expr)* evaluates 'expr' and prints it in LISP-format on the current output.

::

    In> FullForm(2*I*b^2)
    (*
        (Complex 0 2 )
        (^ b 2 ))
    Out> Complex(0,2)*b^2

*Echo(item, item, item, ...)* is a high-level printing routine.

::

    In> Echo({"The square root of four is ", 4*4})
    The square root of four is 16
    Out> True

*PrettyForm(expr)* renders an expression in a nicer way, using ascii art.

::

    In> Taylor(x,0,8) Sin(x)
    Out> x-x^3/6+x^5/120-x^7/5040
    In> PrettyForm(%)

         3    5      7
        x    x      x
    x - -- + --- - ----
        6    120   5040

    Out> True

You can set pretty printing in another way (similar to the init_printing()* method):

::

    In> PrettyPrinter'Set("PrettyForm")

    True

    In> Taylor(x,0,8) Cos(x)

         2    4    6      8
        x    x    x      x
    1 - -- + -- - --- + -----
        2    24   720   40320

*EvalFormula(expr)* displays an evaluation in a nice way, using PrettyPrinter'Set to show 'input=output'.

::

    In> EvalFormula(Taylor(x, 0, 7) Cos(x^2))

                                           4
          /                /  2 \ \       x
    Taylor\ x , 0 , 7 , Cos\ x  / / = 1 - --
                                          2

    Out> True

*TeXForm(expr)* returns a string containing aa LaTeX representation of the Yacas expression expr.

::

    In> TeXForm(Sin(a5)+3*Cos(b5))
    Out> "$\sin a_{5} + 3 \cos b_{5}$"

''''''''''
Plotting
''''''''''

**SymPy**

Pyglet is required to use the plotting function of SymPy in 2d and 3d. Here is an example:

::

    >>> from sympy import symbols, Plot, cos, sin
    >>> x, y = symbols('x y')
    >>> Plot(sin(x*10)*cos(y*5) - x*y)
    [0]: -x*y + sin(10*x)*cos(5*y), 'mode=cartesian'

.. image:: plot.png
   :alt: Output of the plotting example

::

    In[1]: Plot(cos(x*y*10))
    Out[1]: [0]: cos(10*x*y), 'mode=cartesian'

.. image:: plot13.png

::

    In [22]: Plot(1*x**2, [], [x], 'mode=cylindrical') # [unbound_theta,0,2*Pi,40], [x,-1,1,20]
    Out[22]: [0]: x**2, 'mode=cylindrical'

.. image:: plot20.png

**Yacas**

There are two methods of plotting in Yacas:

- Plot2D(f(x), a:b) function
- Plot3DS(f(x,y), a:b, c:d) function

*Plot2D(f, range, option=value)* performs adaptive plotting of one or several functions of one variable in the specified range.

::

    In> Plot2D(Sin(x^4)+Cos(x^4), 0:2*Pi)
    Out> True

.. image:: ya4t.png

*Plot3DS(f, range, option=value)* performs adaptive plotting of two variables in the specified range.


::

    In> Plot3DS(Sin(x*10)*Cos(y*5)-x*y, -1:1, -1:1)
    Out> True

.. image:: ya1t.png

::

    In> Plot3DS(Cos(x*y*10), -1:1, -1:1)
    Out> True

.. image:: ya2t.png

::

    In> Plot3DS(ArcTan(y^3/4-x^2))
    Out> True

.. image:: ya3t.png

''''''''''''
Conclusion
''''''''''''

SymPy aims to be a lightweight normal Python module so as to become a nice open source alternative to Maple/Mathematica. Its goal is to be reasonably fast, easily extended with your own ideas, be callable from Python and could be used in real world problems.
SymPy is perfectly multiplatform, it's small and easy to install and use, since it is written in pure Python (and doesn't need anything else).

Yacas's goal is to be a small system that allows to easily prototype and research symbolic mathematics algorithms. A secondary future goal is to evolve into a complete general purpose CAS.

You can choose to use either SymPy or Yacas, depending on what your needs are. For more information you can go to the official sites of SymPy_ and Yacas_.

.. _SymPy: http://sympy.org/
.. _Yacas: http://yacas.sourceforge.net/homepage.html
