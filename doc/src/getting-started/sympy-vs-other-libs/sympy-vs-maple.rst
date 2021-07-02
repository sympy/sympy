=================
SymPy vs. Maple
=================

|SymPy| vs. |Maple|

.. |SymPy| image:: SymPy.png
.. |Maple| image:: Maple.png

SymPy_ and Maple_ are *Computer algebra systems*.

.. _SymPy: http://sympy.org/
.. _Maple: http://www.maplesoft.com/products/maple/

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
Maple
+++++++

**Maple** is a general-purpose commercial computer algebra system.

Maple was created by the University of Waterloo, Canada and the development began in 1980. Since 1988, Maple has been developed and sold commercially by Waterloo Maple (Maplesoft). The first public release was in 1984.
The latest stable release of Maple is Maple 15.01 (21 June 2011).

Maple is based on a small kernel, written in C, which provides the Maple language. Many libraries from a variety of sources provide functionality to Maple: NAG Numerical Libraries, ATLAS libraries, GMP libraries. Most of the libraries are written in the Maple language and have viewable source code.

Maple is proprietary software restricted by copyright law. You can get Maple for a while if you complete a registration form for an evaluation on their web site. The next step would be to talk to a Maple representative for further information. Maple isn't provided automatically after filling out that form.
If you want to use Maple for more than the evaluation days, then you can choose one of these versions:

* Commercial - $2,275
* Government - $2,155
* Academic - $1245
* Personal Edition - $239
* Student - $99
* Student (12-Month term) - $79

\+ \+: full scientific stack, very functional, fast.

\- \-: very large, expensive, not a library, not open source (proprietary).

++++++++++++++++++++
Sympy (!)= Maple
++++++++++++++++++++

*SymPy* is a cost free open source CAS written entirely in Python released under a modified BSD license while *Maple* is proprietary software released under a proprietary commercial license.

One of the differences between SymPy and Maple is the fact that Maple comes with both a GUI and a command line interface. However, SymPy can use plotting as well, by installing Pyglet.

------------------------
Operating System Support
------------------------

+---------+---------+----------+-------+-----+---------+-----------------------------------+
| System  | Windows | Mac OS X | Linux | BSD | Solaris |               Other               |
+---------+---------+----------+-------+-----+---------+-----------------------------------+
|  SymPy  |   Yes   |    Yes   |  Yes  | Yes |   Yes   |  Any system that supports Python  |
+---------+---------+----------+-------+-----+---------+-----------------------------------+
|  Maple  |   Yes   |    Yes   |  Yes  | No  |   Yes   |                No                 |
+---------+---------+----------+-------+-----+---------+-----------------------------------+

------------------------
Download & Installation
------------------------

Sympy is distributed in various forms. It is possible to download source tarballs and packages from the Google Code page but it is also possible to clone the main Git repository or browse the code online. The only prerequisite is Python since Sympy is Python-based library. It is recommended to install IPython as well, for a better experience.

Maple can be downloaded from the official site after you either get an evaluation version or purchase it.

--------------
Functionality
--------------

+----------+----------+------------+-----------------------------------+---------------------------------------------------------------------------+
|          | Formula  | Arbitrary  |             Calculus              |                                            Solvers                        |
|  System  |          |            +-------------+---------------------+-----------+--------------+-----------------------+------------------------+
|          | editor   | precision  | Integration |Integral transforms* | Equations | Inequalities | Diophantine equations | Differential equations |
+----------+----------+------------+-------------+---------------------+-----------+--------------+-----------------------+------------------------+
|  SymPy   |    No    |    Yes     |    Yes      |        No           |   Yes     |     Yes      |          No           |           Yes          |
+----------+----------+------------+-------------+---------------------+-----------+--------------+-----------------------+------------------------+
|  Maple   |    Yes   |    Yes     |    Yes      |        Yes          |   Yes     |     Yes      |          No           |           Yes          |
+----------+----------+------------+-------------+---------------------+-----------+--------------+-----------------------+------------------------+

+----------+-----------------------+---------+---------+--------------+----------+---------+
|          |        Solvers        | Graph   | Number  | Quantifier   | Boolean  |         |
|  System  +-----------------------+         |         |              |          | Tensors |
|          | Recurrence relations  | theory  | theory  | elimination  | algebra  |         |
+----------+-----------------------+---------+---------+--------------+----------+---------+
|  SymPy   |          Yes          |   No    |   Yes   |     No       |   Yes    |   Yes   |
+----------+-----------------------+---------+---------+--------------+----------+---------+
|  Maple   |          Yes          |   Yes   |   Yes   |     No       |   No     |   Yes   |
+----------+-----------------------+---------+---------+--------------+----------+---------+

\* Will be available in SymPy 0.7.2

**Note**

This document contains some examples from Maple's documentation and are under a different license than SymPy.

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

However, in Maple, both ^ and \*\* mean exponentiation:

::

    > (x+1)^2;
           2
    (x + 1)
    > (x+1)**2;
           2
    (x + 1)

Another difference between SymPy and Maple is that you have to define symbols in SymPy before you can use them, while in Maple it isn't necessary.

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

**Maple**

::

    > x**2 + 2*x + 1;
    x^2 + 2*x + 1

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

**Maple**

Here are some examples of algebra in Maple:

*Expand(expr, expr1, expr2,..., exprn)* expands an expression. This method is equal to the *apart(expr, x)* method from SymPy.

::

    > expand((x^2+x+1)/(x+2));
      2
     x        x       1
    ----- + ----- + -----
    x + 2   x + 2   x + 2

The *evalf()* method evaluates expressions:

::

    > evalf(Pi)
    3.141592654

    > evalf(sqrt(2)*Pi)
    4.442882938

*Divide(a, b, 'q', options)* does a check for exact divisibility for polynomials with algebraic number coefficients.

::

    > with(Algebraic):
    > Divide(x^2+1, x+I, 'q');
    true

    > q;
    x-I

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

**Maple**

The *limit(f, x=a, dir)* function attempts to compute the limiting value of f as x approaches a.

::

    > limit(sin(x)/x, x = 0);
    1

Limit of f(x)= 2*x+1 as x -> 5/2

::

    > limit(2*x+1, x = 5/2);
    6

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

**Maple**

Here are some examples of differentiation:

::

    > diff(cos(x^3), x)

.. image:: o17.png

::

    > diff(1/tan(x), x)

.. image:: o18.png

::

    > diff(1/cot(x^2), x)

.. image:: o19.png

*BesselJ(v, x)* is the Bessel function of the first kind. It satisfies Bessel's equation: x^2*y''+x*y'+(x^2-v^2)*y = 0

::

    > BesselJ(v, x)
    BesselJ(v, x)

    > diff(BesselJ(v, x), x)

.. image:: o20.png

::

    > series(BesselJ(5, x), x)

.. image:: o21.png

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

**Maple**

The *taylor(expr, x=a, n)* command computes the order n Taylor series expansion of expr, with respect to the variable x, about the point a.

::

    > taylor(cos(x), x = 0, 14)

.. image:: o15.png

::

    > taylor(1/cos(x^2), x = 0, 14)

.. image:: o16.png

*order(expr)* determines the truncation order of a series:

::

    > order(taylor(1/cos(x^2), x = 0, 14))+28
    42

    > order(taylor(1/cos(x^2), x = 0, 14))+x^2
    14+x^2

    > series(1/(1-x), x)
    1+x+x^2+x^3+x^4+x^5+O(x^6)		(1)

    > order((1))
    6

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

**Maple**

The *int(expr, x)* calling sequence computes an indefinite integral of the expr with respect to the variable x.

* polynomial functions:

::

    > f := x^2+2*x+4:
    > int(f, x)

.. image:: o1.png

* rational functions:

::

    > f := (x+1)/(x^2+4*x+4):
    > int(f, x)

.. image:: o2.png

* exponential-polynomial functions:

::

    > f := 5*x^2*exp(x)*sin(x):
    > int(f, x)

.. image:: o3.png

* non-exponential integrals:

::

    > f := exp(-x^2)*erf(x):
    > int(f, x)

.. image:: o4.png

The *int(expr, x=a..b)* calling sequence computes the definite integral of the expr with respect to the variable x on the interval from a to b.

::

    > int(x^2*cos(x), x = 0 .. (1/2)*Pi)

.. image:: o5.png

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

**Maple**

::

    > exp(I*2*x);

.. image:: o12.png

To return a complex number in 'a+bI' format in Maple, there is the *evalc(expr)* command, which is a symbolic evaluator over the complex field.

::

    > evalc(exp(2*I*x));

.. image:: cn1.png

::

    > assume(x::real);
    > evalc(exp((2*I)*x));

.. image:: cn2.png

::

    > assume(x::real);
    > evalc(exp(-2+(3*I)*x));

.. image:: cn3.png

Complex numbers division in Maple:

::

    > (2+3*I)/(3+7*I);

.. image:: cn4.png

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

    In [13]: gamma(b+2).series(b, 0, 3)
    Out[13]:
                            2  2             2  2
                           π ⋅b    EulerGamma ⋅b                2    ⎛ 3⎞
    1 + b - EulerGamma⋅b + ───── + ────────────── - EulerGamma⋅b  + O⎝b ⎠
                             12          2

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

**Maple**

**trigonometric**

::

    > expand(cos(x-y))
    cos(x)*cos(y)+sin(x)*sin(y)

    > expand(cos(2*x))
    2*cos(x)^2-1

    > sinh(I*x^2)
    I*sin(x^2)

    > simplify(sin(x)^2*cos(y)^2+cos(x)^2*sin(y)^2, trig)

.. image:: o11.png

**zeta function**

::

    > Zeta(5, 2)

.. image:: o9.png

::

    > Zeta(1, 1/2)

.. image:: o10.png

::

    > Zeta(2.2)
    1.490543257

**factorials and gamma function**

::

    > a!
    a!

    > 10!
    3628800

    > GAMMA(1.0+2.5*I, 2.0+3.5*I)
    0.1314614269e-1+0.6253182683e-2*I

**polynomials**

::

    > ChebyshevT(8, x)
    ChebyshevT(8, x)

    > series((4),'ChebyshevT')

.. image:: o6.png

::

    > LegendreP(3, x)
    LegendreP(3, x)

    > series((6),'LegendreP')

.. image:: o7.png

::

    > HermiteH(3, x^2)
    HermiteH(3, x^2)		(8)

    > series((8),'HermiteH')

.. image:: o8.png

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

**Maple**

The *D* and *diff* commands can both compute derivatives. The *D* operator computes derivatives of operators, while *diff* computes derivatives of expressions.

::

    > D(ln)(x)
    1/x		(4)

    > convert((4), diff)
    1/x

    > diff(sin(x), x)
    cos(x)

The *dsolve(ODE, y(x), options)* command solves ordinary differential equations (ODEs):

::

    > dsolve(diff(f(x), x)+f(x), f(x))

.. image:: o26.png

In this example, we define a derivative using the diff command and solve the ODE.

::

    > ode := diff(y(x), x, x) = 2*y(x)+1

.. image:: o27.png

::

    > dsolve(ode)

.. image:: o28.png

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

**Maple**

::

    > solve(x^3+2*x^3-1, x)

.. image:: o13.png

::

    > solve({x+2*y = 3, y+1/x = 1}, [x, y])

.. image:: o14.png

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

**Maple**

You can create a matrix in Maple by using the Matrix palette or by using Maple's matrix notation.
Here are a few examples:

::

    > <<1, 0, 0>|<0, 1, 0>|<0, 0, 1>>
    (1)

.. image:: ma1.png

To assign a name to a matrix you have to type *name := * and use an equation label to refer to the matrix by typing Ctrl+L. In the *Insert Label* box type 1.

::

    > A := (1)
    (2)

.. image:: ma1.png

To perform matrix calculations, you must use the context menu.
Right click on the matrix and select *Standard Operations>Determinant* to find the determinant.

::

    > <<2, 4, 1>|<5, 7, 0>|<6, 10, 3>>

.. image:: ma2.png

To compute the transpose of a matrix you have to go to *Standard Operations>Transpose*.

::

    > <<1, 3>|<2+I, 4>>

.. image:: ma3.png

The function *hadamard(A)* computes a bound on the maxnorm of det(A) where A is an n by n matrix.

::

    > with(LinearAlgebra):
    > A := matrix(2, 2, [1, 3, 18, 3]);

.. image:: o22.png

::

    > B := matrix(2, 2, [10, 4, 20, 6]);

.. image:: o23.png

::

    > hadamard(A);

.. image:: o24.png

::

    > hadamard(B);

.. image:: o25.png

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

**Maple**

The commands in the geometry module enable you to work in two-dimensional Euclidean geometry.

Define a point (*point(P,Px,Py)*):

::

    > with(geometry):
    > point(A, 2, 3);
    A

The *detail(P)* command returns a detailed expression of the point P.

::

    > detail(A);

.. image:: o29.png

Check if points are collinear:

Note that it is necessary to define at least three points for the *AreCollinear()* command because it has only the following syntax: *AreCollinear(P, Q, R, cond)*, where P, Q, R are three previously defined points and cond is an optional name.

::

    > with(geometry):
    > point(A, 0, 0), point(B, 3, 1), point(C, 5, 5), point(F, 2, 2);
    > AreCollinear(A, B, C);
    false
    > AreCollinear(A, C, F);
    true

Define a segment (*segment(seg, [P1, P2])*):

::

    > with(geometry):
    > point(A, 1, 2), point(B, 2, -3);
    A, B
    > segment(AB, [A,B])
    AB
    > DefinedAs(AB)
    [A,B]		(3)
    > map(coord
    [ [1, 2], [2, -3] ]

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

**Maple**

The *match(expr = pattern, v, 's')* calling sequence returns true if it can match expr to pattern for some values of the variables (excluding the main variable, v). Otherwise, it returns false.

::

    > match(10*x^3 = y*x^3, x, 's')
    true

    > s;
    {y = 10}

    > match(x^4 = y*x^q, x, 's')
    true

    > s;
    {q = 4, y = 1}

''''''''''
Printing
''''''''''

**SymPy**

There are many ways of printing mathematical expressions.
Two of the most common methods are:

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

**Maple**

There are several methods to print expressions in Maple:

* print(e1, e2, ...) command
* lprint(expr1, expr2,...) command
* fprintf, sprintf, nprintf, printf commands
* 2d ascii pretty printing in command line

The *print(e1, e2, ...)* command is equal to the *pprint()* method from SymPy.

::

    > print(c^3);

.. image:: o30.png

::

    > print(2/c);

.. image:: o31.png

::

    > print('int(c^2+2*c, c)');

.. image:: o32.png

The *lprint(expr1, expr2, ...)* performs linear printing of expressions and is the same as the standard printing (the *print()* method) from SymPy.

::

    > lprint(sin(x)^2+cos(x)^2);

.. image:: o33.png

::

    > lprint(int(sin(x+y)/(x-y), x))

.. image:: o34.png

Note that in command-line Maple, the output of *lprint* cna be cut and paste into a Maple session while pretty-printed output cannot.

There are other printing commands that have similar syntax to the print command:

* fprintf - prints expressions to a file or pipe based on a format string
* sprintf - prints expressions to a string based on a format string
* nprintf - prints expressions to a name based on a format string
* printf - prints expressions to a default stream based on a format string

The *fprintf()* command returns a count of the number of characters written.

::

    > fd := fopen("temp_file", WRITE);
    1
    > fprintf(fd, "x = %d, y = %g", 2, 1.5);
    14
    > fclose(fd);

The *nprintf()* command is the same as *sprintf()*, except that it returns a Maple symbol (a simple name) instead of a string.

::

    > sprintf("%o %x", 2805, 2805); or nprintf("%o %x", 2805, 2805);

.. image:: o35.png
.. image:: o36.png

::

    > printf("%-2.5s:%2.5s:%2.5s", S, Sym, SymPyV);

.. image:: o37.png

Maple also has a 2d ASCII pretty printer that is used with the command line interface.

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

**Maple**

There are two main methods to generate plots in Maple:

* plot(f, x=x0..x1) command to create a 2d plot
* plot3d(expr, x=a..b, y=c..d, opts) command for 3d plotting

Here are some examples:

::

    > plot(BesselJ(1, x), x = -30 .. 30);

.. image:: plot15.png

::

    > plot3d({x+2*y, sin(x*y)}, x = -Pi .. Pi, y = -Pi .. Pi);

.. image:: plot16.png

::

    > plot3d(cos(10*x*y), x = -1 .. 1, y = -1 .. 1);

.. image:: plot17.png

::

    > plot3d(height, angle = 0 .. 2*Pi, height = -5 .. 5, coords = cylindrical, title = CONE);

.. image:: plot18.png

::

    > plot3d(r*cos(theta), r = 0 .. 10, theta = 0 .. 2*Pi, coords = cylindrical, orientation = [100, 71], axes = NONE);

.. image:: plot19.png

''''''''''''
Conclusion
''''''''''''

SymPy aims to be a lightweight normal Python module so as to become a nice open source alternative to Maple. Its goal is to be reasonably fast, easily extended with your own ideas, be callable from Python and could be used in real world problems.
Another advantage of SymPy compared to Maple is that since it is written in pure Python (and doesn't need anything else), it is perfectly multiplatform, it's small and easy to install and use.

You can choose to use either SymPy or Maple, depending on what your needs are. For more information you can go to the official sites of SymPy_ and Maple_.

.. _SymPy: http://sympy.org/
.. _Maple: http://www.maplesoft.com/products/maple/