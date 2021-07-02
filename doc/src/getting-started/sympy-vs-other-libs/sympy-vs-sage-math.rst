====================
SymPy vs. SageMath
====================

|SymPy| vs. |SageMath|

.. |SymPy| image:: SymPy.png
.. |SageMath| image:: Sage.png

SymPy_ and SageMath_ are *Computer algebra systems*.

.. _SymPy: http://sympy.org/
.. _SageMath: http://www.sagemath.org/

**Computer Algebra System**
    A software program that facilitates symbolic mathematics.
    The core functionality of a CAS is manipulation of mathematical expressions in symbolic form.

+++++++
SymPy
+++++++

**SymPy** is a Python library for symbolic computation that aims to become a full-featured computer algebra system and to keep the code simple to promote extensibility and comprehensibility.

SymPy was started by Ondřej Čertík in 2005 and he wrote some code in 2006 as well. In 11 March 2007, SymPy was realeased to the public.
The latest stable release of SymPy as of time of writing is 1.1.1. There are more than 630 contributors to SymPy.

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

++++++
SageMath
++++++

**SageMath** (nickname **Sage**) is a mathematical software with features covering many aspects of mathematics, including algebra, number theory and calculus.

Sage was created by William A. Stein and the development began in 2005. The first public release was in the same year (24 February 2005).
The latest stable release of Sage is Sage 8.7 (23 March 2019). There are more than 150 contributors to Sage.

Sage can be used in several ways:

- Cloud (in the browser)
- Notebook graphical interface
- Interactive command line
- Programs (writing interpreted programs)
- Scripts (writing stand-alone Python scripts that use the Sage library)

The overall goal of Sage is to create a viable, free, open-source alternative to Maple, Mathematica, Magma and MATLAB.
Sage integrates many specialized mathematics software into a common interface, for which a user needs to know only Python.

However, Sage is not pure Python, since it makes use of many external packages written in many other languages and it is actually a collection of dozens different computer algebra systems.

\+ \+: full scientific stack, very functional, fast, large community.

\- \-: very large, not a library, complicated design.

++++++++++++++++++++
Sympy (!)= SageMath
++++++++++++++++++++

Both *SymPy* and *Sage* are cost free open source CASes written entirely in Python (although Sage contains Cython code as well), the first one being released under a modified BSD license and the latter under a GNU GPL license.

One of the differences between SymPy and Sage is the fact that Sage comes with a GUI - the Jupyter Notebook. This Notebook is useful because it allows you to write and run code, display 2d and 3d plots and organize and share your work. The Notebook is run by typing *./sage -n jupyter* after you cd into the Sage directory.
However, SymPy can use plotting as well, by installing Pyglet or matplotlib.

------------------------
Operating System Support
------------------------

+---------+---------+----------+-------+-----+---------+-------------------------------------+
| System  | Windows | Mac OS X | Linux | BSD | Solaris |                Other                |
+---------+---------+----------+-------+-----+---------+-------------------------------------+
|  SymPy  |   Yes   |    Yes   |  Yes  | Yes |   Yes   |  Any system that supports Python    |
+---------+---------+----------+-------+-----+---------+-------------------------------------+
|SageMath |   Yes   |    Yes   |  Yes  | No  |   Yes   |                                     |
+---------+---------+----------+-------+-----+---------+-------------------------------------+

------------------------
Download & Installation
------------------------

Sympy is distributed in various forms. It is possible to download source tarballs and packages from the Google Code page but it is also possible to clone the main Git repository or browse the code online. The only prerequisite is Python since Sympy is Python-based library. It is recommended to install IPython as well, for a better experience.

Sage is completely usable in the browser. Sage can also be installed either from a pre-built binary tarball or from source. The binary method is fastest and has fewest prerequisites. The source method may take a long time (from 1 hour to 14 days), depending on the computer (it took about two weeks to build Sage on the T-Mobile G1 Android telephone).

--------------
Functionality
--------------

+----------+----------+------------+-----------------------------------+---------------------------------------------------------------------------+
|          | Formula  | Arbitrary  |             Calculus              |                                            Solvers                        |
|  System  |          |            +-------------+---------------------+-----------+--------------+-----------------------+------------------------+
|          | editor   | precision  | Integration |Integral transforms  | Equations | Inequalities | Diophantine equations | Differential equations |
+----------+----------+------------+-------------+---------------------+-----------+--------------+-----------------------+------------------------+
|  SymPy   |    No    |    Yes     |    Yes      |        Yes          |   Yes     |     Yes      |          Yes          |           Yes          |
+----------+----------+------------+-------------+---------------------+-----------+--------------+-----------------------+------------------------+
| SageMath |    Yes   |    Yes     |    Yes      |        Yes          |   Yes     |     Yes      |          No           |           Yes          |
+----------+----------+------------+-------------+---------------------+-----------+--------------+-----------------------+------------------------+

+----------+-----------------------+---------+---------+--------------+----------+---------+
|          |        Solvers        | Graph   | Number  | Quantifier   | Boolean  |         |
|  System  +-----------------------+         |         |              |          | Tensors |
|          | Recurrence relations  | theory* | theory  | elimination  | algebra  |         |
+----------+-----------------------+---------+---------+--------------+----------+---------+
|  SymPy   |          Yes          |   No    |   Yes   |     No       |   Yes    |   Yes   |
+----------+-----------------------+---------+---------+--------------+----------+---------+
| SageMath |          Yes          |   Yes   |   Yes   |     Yes      |   Yes    |   Yes   |
+----------+-----------------------+---------+---------+--------------+----------+---------+

\* It is recommended that people use NetworkX_ for Graph Theory.

.. _NetworkX: https://networkx.github.io

'''''''''''''''''''''''''
Some syntax differences
'''''''''''''''''''''''''

*Sage* and *SymPy* may look very similar, but those are two very different systems with completely different internal design, non-overlapping features sets (e.g. Sage is very good at number theory and abstract algebra, but SymPy has sophisticated pretty printing and code generation frameworks) and quite different semantics.

SymPy uses Python constructs only. Here is an example:

::

    >>> 2/7        # Python evaluates this to 0
    0

    >>> from __future__ import division         # We obtain a different result if we import division from __future__
    >>> 2/7
    0.285714285714

In Sage, the example returns a Rational:

::

    sage: 2/7
    2/7

To obtain a Rational in SymPy, one of these methods must be used:

::

    >>> from sympy import Rational
    >>> Rational(2, 7)
    2/7

    >>> from sympy import S
    >>> S(2)/7
    2/7

In SymPy, to raise something to a power, you must use \*\*, not ^ as the latter uses the Python meaning, which is xor.

However, in Sage, both \*\* and ^ can be used to perform exponentation. This shows that Sage has a modified version of Python.

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

    sage: (x+1)^2
    (x+1)^2

    sage: (x+1)**2
    (x+1)^2

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

**Sage**

Here are some examples of algebra in Sage:

::

    sage: 1/( (x**2+2*x+1)*(x**2-1) )
    1/((x^2 - 1)*(x^2 + 2*x + 1))

    sage: expand((x-1)^2)
    x^2 - 2*x + 1

    sage: f = I + x - x
    sage: simplify(f)
    I

    sage: f = (cos(x)*sin(y))/sin(y)+(sin(x)*cos(y))/sin(x)
    sage: simplify(f)
    cos(x) + cos(y)

In Sage, to return the exact value of expressions, *n()*, *.n(digits)* and *numerical_approx(var, prec)* are used:

::

    sage: n(pi)
    3.14159265358979

    sage: N(sqrt(2)*pi, digits=50)
    4.4428829381583662470158809900606936986146216893757

    sage: numerical_approx(pi, prec=200)
    3.1415926535897932384626433832795028841971693993751058209749

    sage: f = x^(-2*x)
    sage: f.integral(x, 1, +Infinity)
    integrate(x^(-2*x), x, 1, +Infinity)
    sage: show(integrate(x^(-2*x), x, 1, +Infinity))

.. image:: int4.png

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

**Sage**

The limit() function has this syntax: *limit(ex, dir=None, taylor=False, algorithm=’maxima’, **argv)*:

::

    sage: x = var('x')
    sage: f = (1+1/x)^x
    sage: f.limit(x = oo)
    e

    sage: f.limit(x = 5)
    7776/3125

    sage: f.limit(x = I, taylor=True)
    (-I + 1)^I

    sage: limit(2*x+1, x=5/2)       # 5/2 is Rational in Sage by default
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

**Sage**

::

    sage: var ('x')        # declaration of variable x
    x
    sage: f = cos(x**3); f
    cos(x^3)
    sage: f.diff(x)
    -3*x^2*sin(x^3)
    sage: show(f)          # show() is one of the pprint functions in Sage

.. image:: so1.png

::

    sage: show(f.diff(x))

.. image:: so2.png

::

    sage: f = atan(2*x); f
    arctan(2*x)
    sage: f.diff(x)
    2/(4*x^2 + 1)
    sage: show(f)

.. image:: so3.png

::

    sage: show(f.diff(x))

.. image:: so4.png

::

    sage: f = 1/tan(x); f
    1/tan(x)
    sage: f.diff(x)
    -(tan(x)^2 + 1)/tan(x)^2
    sage: show(f)

.. image:: so5.png

::

    sage: show(f.diff(x))

.. image:: so6.png

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

**Sage**

Series expansion can be used with the *taylor(f, *args)* function:

::

    sage: taylor(cos(x), x, 0, 14)
    -1/87178291200*x^14 + 1/479001600*x^12 - 1/3628800*x^10 + 1/40320*x^8 - 1/720*x^6 + 1/24*x^4 - 1/2*x^2 + 1

    sage: taylor(1/cos(x**2), x, 0, 14)
    61/720*x^12 + 5/24*x^8 + 1/2*x^4 + 1

    sage: var('x, k')
    (x, k)
    sage: taylor(sqrt (1 - k^2*sin(x)^2), x, 0, 6)
    -1/720*(45*k^6 - 60*k^4 + 16*k^2)*x^6 - 1/24*(3*k^4 - 4*k^2)*x^4 - 1/2*k^2*x^2 + 1

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

**Sage**

Sage can integrate some simple functions on its own:

- polynomial functions:

::

    sage: f = x^2 + 2*x + 4     # ** and ^ represent the same thing
    sage: f.integral(x)
    1/3*x^3 + x^2 + 4*x

- rational functions:

::

    sage: f = (x+1)/(x^2 + 4*x + 4)
    sage: f.integral(x)
    1/(x + 2) + log(x + 2)

- exponential-polynomial functions:

::

    sage: f = 5*x^2 * exp(x) * sin(x)
    sage: f.integral(x)
    5/2*(x^2 - 1)*e^x*sin(x) - 5/2*(x^2 - 2*x + 1)*e^x*cos(x)

- non-elementary integrals:

::

    sage: f = exp(-x**2)*erf(x)
    sage: f.integral(x)
    1/4*sqrt(pi)*erf(x)^2

The output of |integral2| in Sage is:

.. |integral2| image:: int3.png

::

    sage: f = x^2 * cos(x)
    sage: f.integral(x, 0, pi/2)
    1/4*pi^2 - 2


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

**Sage**

This is an example of a complex number in Sage:

::

    sage: C = ComplexField()
    sage: I = C.0
    sage: b = 15/10 + 25/10*I
    sage: b
    1.50000000000000 + 2.50000000000000*I

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

**Sage**

**trigonometric**

::

    sage: (cos(90)+sin(30)).n(100)
    -1.4361052402220319423532262217

    sage: (cos(3) + tan(3)).n(5)
    -1.1

    sage: sinh(arccosh(x), hold=True).simplify()
    sqrt(x - 1)*sqrt(x + 1)

**zeta function**

The example below returns the Riemann zeta function evaluated at a complex number:

::

    sage: i = ComplexField(30).gen()
    sage: z = 1 + i
    sage: z.zeta()
    0.58215805981 - 0.92684856430*I
    sage: zeta(z)
    0.58215805981 - 0.92684856430*I

**factorials and gamma function**

The *factorial(*args, coerce=True, hold=False, dont_call_method_on_arg=False)* function return the factorial of n (the output is an integer or a symbolic expression).

::

    sage: x = var('x')
    sage: factorial(10)
    3628800

The *gamma()* function is used for other nonnegative numbers that are not integers:

::

    sage: x = var('x')
    sage: factorial(3/4)
    gamma(7/4)
    sage: factorial(2.3)
    2.68343738195577

These examples return the Gamma function and the incomplete form of it, evaluated for a complex number:

::

    sage: i = ComplexField(30).0
    sage: (1+i).gamma()
    0.49801566824 - 0.15494982828*I

    sage: C, i = ComplexField(30).objgen()
    sage: (1+i).gamma_inc(2 + 3*i)
    0.0020969148645 - 0.059981913655*I

**polynomials**

::

    sage: x = polygen(QQ, 'x')
    sage: chebyshev_U(2, x)
    4*x^2 - 1

Special functions like chebyshev or bessel have only numerical use in Sage. For symbolic use, the Maxima interface included in Sage must be used directly:

::

    sage: maxima.eval("f:bessel_y(v, w)")
    'bessel_y(v,w)'

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

**Sage**

You can use Sage to investigate ordinary differential equations. To solve the equation x'+x-1=0:

::

    sage: t = var('t')
    sage: x = function('x', t)
    sage: d = diff(x, t) + x - 1
    sage: desolve(d, [x, t])
    (c + e^t)*e^(-t)

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

**Sage**

The *solve()* function solves equations. To use it, you must specify some variables, then the argumens to solve are an equation (or a system of equations), together with the variables for which to solve. The function is the same as the solve() from SymPy:

::

    sage: x = var('x')
    sage: solve(x^3 + 2*x^3 - 1, x)
    [x == 1/6*(I*sqrt(3) - 1)*3^(2/3), x == 1/6*(-I*sqrt(3) - 1)*3^(2/3), x == 1/3*3^(2/3)]

    sage: x = var('x')
    sage: solve( (x^2 + 4*y^2 -2, -10*x + 2*y -15), (x, y) )
    [ [x == -1/101*I*sqrt(23) - 150/101, y == -5/101*I*sqrt(23) + 15/202], [x == 1/101*I*sqrt(23) - 150/101, y == 5/101*I*sqrt(23) + 15/202] ]

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

**Sage**

Sage provides standard constructions from linear algebra:

::

    sage: A = Matrix([ [1, 0 , 0], [0, 1, 0], [0, 0, 1] ])
    [1 0 0]
    [0 1 0]
    [0 0 1]

*det(x)* returns the determinant of x:

::

    sage: M = MatrixSpace(QQ, 3, 3)
    sage: A = M([2, 5, 6, 4, 7, 10, 1, 0, 3])
    sage: det(A)
    -10

It is not possible to define a matrix with various types of numbers.

::

    sage: M = MatrixSpace(CC, 2, 2)   # Complex numbers matrix
    sage: A = M( [1,2+I,3,4] )      # both integer and complex numbers
    sage: transpose(A)
    [1.00000000000000                           3.00000000000000]
    [2.00000000000000 + 1.00000000000000*I      4.00000000000000]

    sage: M = MatrixSpace(ZZ, 2, 2)   # integer numbers
    sage: A = M( [1,2,3,4] )
    sage: transpose(A)
    [1 3]
    [2 4]

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

**Sage**

The *point(points, **kwds)* function return either a 2d or 3d point or sum of points. Here are some examples:

::

    sage: point([(1,2), (1,3), (2,2)])

.. image:: points1.png
   :alt: Points of a triangle

::

    sage: point([(cos(theta), sin(theta)) for theta in srange(0, 2*pi, pi/8)]).show(frame=True)

.. image:: points2.png
   :alt: The trigonometric circle highlighted by points

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

**Sage**

The *wild(n=0)* function returns the n-th wild-card for pattern matching and substitution:

::

    sage: x,y = var('x,y')
    sage: w0 = SR.wild(0); w1 = SR.wild(1)
    sage: pattern = sin(x)*w0*w1^2; pattern
    $0*$1^2*sin(x)
    sage: f = atan(sin(x)*3*x^2); f
    arctan(3*x^2*sin(x))
    sage: f.has(pattern)
    True
    sage: f.subs(pattern == x^2)
    arctan(x^2)

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

*Pretty printing* is a nice ascii-art printing with the help of a *pprint* function:

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

**Sage**

There are three common ways to print expressions in Sage:

- Standard printing
- The *show()* function
- Pretty printing

*Standard printing* in Sage is almost identical with the printing from SymPy (Note that there is no print() function in Sage):

::

    sage: c = var('c')
    c^3
    sage: 2/c
    2/c
    sage: integrate(c**2+2*c, c)
    1/3*c^3 + c^2

The *show()* function is the same as the *pprint()* function from SymPy:

::

    sage: f = x^(-2*x)
    sage: f.integral(x, 1, +Infinity)
    integrate(x^(-2*x), x, 1, +Infinity)
    sage: show(integrate(x^(-2*x), x, 1, +Infinity))

.. image:: int4.png

The *pretty_print()* function can be considered similar to the *show()* function, but the output is slightly different (pretty_print() displays the integral on one line, whereas show() displays it on three lines.):

::

    sage: pretty_print(integrate(x^(-2*x), x, 1, +Infinity))

.. image:: int5.png

''''''''''
Plotting
''''''''''

**SymPy**

Pyglet or matplotlib are required to use the plotting function of SymPy in 2d and 3d. Here is an example using pyglet:

::

    >>> from sympy import symbols, Plot, cos, sin
    >>> x, y = symbols('x y')
    >>> Plot(sin(x*10)*cos(y*5) - x*y)
    [0]: -x*y + sin(10*x)*cos(5*y), 'mode=cartesian'

.. image:: plot.png
   :alt: Sympy output of the 3-d plotting example using pyglet

And the same example using matplotlib (though with finer sampling):

::

    >>> from sympy import symbols, cos, sin
    >>> from sympy.plotting import plot3d
    >>> x, y = symbols('x y')
    >>> plot3d(sin(x*10)*cos(y*5) - x*y, (x, -1, 1), (y, -1, 1),
               nb_of_points_x=100, nb_of_points_y=100)

.. image:: https://raw.githubusercontent.com/wiki/MOBle/sympy/plot3d_matplotlib.png
   :alt: Sympy output of the 3-d plotting example using matplotlib

**Sage**

In Sage, you can produce filled-in shapes by creating a list of points (L in the example below) and then use the *polygon* command to plot the shape with boundary formed by those points.
For example, here is a blue hypotrochoid. By typing *show(p, axes=false)*, you can see this without any axes. Note that it is possible to add text to a plot.

::

    sage: L = [[6*cos(pi*i/100)+5*cos((6/2)*pi*i/100), 6*sin(pi*i/100)-5*sin((6/2)*pi*i/100)] for i in range(200)]
    sage: p = polygon(L, rgbcolor=(1/8,1/4,1/2))
    sage: t = text("hypotrochoid", (5,4), rgbcolor=(1,0,0))
    sage: show(p+t, axes=false)

.. image:: plot2.png
   :alt: Hypotrochoid (axes are disabled)

The function sin(x)/x has the following plot:

::

    sage: plot(sin(x)/x, x, -100, 100).show(ymin=-1)

.. image:: plot3.png
   :alt: Output of sin(x)/x

''''''''''''
Conclusion
''''''''''''

SymPy and Sage are trying to become nice open source alternatives to Maple/Mathematica. Their goal is to be reasonably fast, easily extended with your own ideas, be callable from Python and could be used in real world problems. SymPy uses a different approach to achieve this goal, because it aims to be a lightweight normal Python module, whereas Sage aims to glue together every useful open source mathematics software package (that is why SymPy is included in Sage by default since version 2.7) and provide a transparent interface to all of them.
Another advantage of SymPy is that since it is written in pure Python (and doesn't need anything else), it is perfectly multiplatform, it's small and easy to install and use.

You can choose to use either SymPy or Sage, depending on what your needs are. For more information you can go to the official sites of SymPy_ and Sage_.

.. _SymPy: http://sympy.org/
.. _Sage: http://www.sagemath.org/