=================
 SymPy vs. Axiom
=================

|SymPy| vs. |Axiom|

.. |SymPy| image:: SymPy.png
.. |Axiom| image:: axiom.png

SymPy_ and Axiom_ are *Computer algebra systems*.

.. _SymPy: http://sympy.org/
.. _Axiom: http://axiom-developer.org/

**Computer Algebra System**
    A software program that facilitates symbolic mathematics.
    The core functionality of a CAS is manipulation of mathematical expressions in symbolic form.

+++++++
SymPy
+++++++

**Sympy** is a Python library for symbolic computation that aims to become a full-featured computer algebra system and to keep the code simple to promote extensibility and comprehensibility.

SymPy was started by Ondřej Čertík in 2005 and he wrote some code in 2006 as well. In 11 March 2007, SymPy was realeased to the public. The latest stable release of SymPy is 1.0 (March 3, 2016). As of beginning of December 2001 there have been over 150 people who contributed at least one commit to SymPy.

SymPy can be used:

- Inside Python, as a library.
- As an interactive command line, using IPython.

SymPy is entirely written in Python and does not require any external libraries, but various programs that can extend its capabilities can be installed:

- gmpy, Cython --> speed improvement.
- Pyglet, Matplotlib --> 2d and 3d plotting.
- IPython --> interactive sessions.

SymPy is available online at `SymPy Live`_. The site was developed specifically for SymPy. It is a simple web shell that looks similar to iSymPy under the standard Python interpreter. SymPy Live uses Google App Engine as computational backend.

.. _`SymPy Live`: http://live.sympy.org/


+++++++
Axiom
+++++++

**Axiom** is a general purpose Computer Algebra system. It is useful for research and development of mathematical algorithms. It consists of an interpreter environment, a compiler and a library, which defines a strongly typed, mathematically (mostly) correct type hierarchy.

Axiom has been in development since 1971, originally as Scratchpad by researchers at IBM under the direction of Richard Dimick Jenks. Other key early developers were Barry Trager, Stephen Watt, James Davenport, Robert Sutor, and Scott Morrison.
In the 1990s it was sold to NAG and given its current name. In 2001 it was withdrawn from the market and re-released under the Modified BSD License. Since then, the project lead developer has been Tim Daly.
In 2007, Axiom was forked into two different open-source projects: OpenAxiom, and FriCAS. Last stable release was in August 2014.

- Axiom is a literate program.
- The source code is becoming available in a set of volumes which are available on the `Axiom Developer`_ website.
These volumes contain the actual source code of the system.

.. _`Axiom Developer`: https://axiom-developer.org/

++++++++++++++++++++++++++++++++++++
Differences between Sympy and Axiom
++++++++++++++++++++++++++++++++++++

-------------------------
Operating System Support
-------------------------

+----------------+---------+----------+-------+-----+---------+-------------------------------------+
| System         | Windows | Mac OS X | Linux | BSD | Solaris |                Other                |
+----------------+---------+----------+-------+-----+---------+-------------------------------------+
|  SymPy         |   Yes   |    Yes   |  Yes  | Yes |   Yes   |  Any system that supports Python    |
+----------------+---------+----------+-------+-----+---------+-------------------------------------+
|  Axiom         |   Yes   |    Yes   |  Yes  | No  |   No    |                  No                 |
+----------------+---------+----------+-------+-----+---------+-------------------------------------+

-------------------------
Download & Installation
-------------------------

Sympy is distributed in various forms. It is possible to download source tarballs and packages from the Google Code page but it is also possible to clone the main Git repository or browse the code online. The only prerequisite is Python since Sympy is Python-based library. It is recommended to install IPython as well, for a better experience.

Axiom can be downloaded and installed from the `Axiom Developer`_ website, also available on Github.

.. _`Axiom Developer`:https://axiom-developer.org/

---------------
Functionality
---------------

As you can see in the tables below, there are some differences between the available features in both mathematical systems.

+---------------+----------+------------+-----------------------------------+---------------------------------------------------------------------------+
|               | Formula  | Arbitrary  |             Calculus              |                                  Solvers                                  |
|  System       |          |            +-------------+---------------------+-----------+--------------+-----------------------+------------------------+
|               | editor   | precision  | Integration |Integral transforms  | Equations | Inequalities | Diophantine equations | Differential equations |
+---------------+----------+------------+-------------+---------------------+-----------+--------------+-----------------------+------------------------+
|  SymPy        |    No    |    Yes     |    Yes      |        No           |   Yes     |     Yes      |          No           |           Yes          |
+---------------+----------+------------+-------------+---------------------+-----------+--------------+-----------------------+------------------------+
|  Axiom        |   No     |    Yes     |    Yes      |       No            |     Yes   |    Yes       |      Yes          |             Yes          |
+---------------+----------+------------+-------------+---------------------+-----------+--------------+-----------------------+------------------------+

+---------------+-----------------------+---------+---------+--------------+----------+---------+
|               |        Solvers        | Graph   | Number  | Quantifier   | Boolean  |         |
|  System       +-----------------------+         |         |              |          | Tensors |
|               | Recurrence relations  | theory  | theory  | elimination  | algebra  |         |
+---------------+-----------------------+---------+---------+--------------+----------+---------+
|  SymPy        |          Yes          |   No    |   Yes   |     No       |   Yes    |   Yes   |
+---------------+-----------------------+---------+---------+--------------+----------+---------+
| Axiom         |          Yes          |   No    |   Yes   |     Yes      |   Yes    |   Yes   |
+---------------+-----------------------+---------+---------+--------------+----------+---------+

-------------------------------------
Some Features and syntax differences
-------------------------------------

One very important aspect of Axiom is that it uses types extensively. Everything you do in Axiom takes an input of one of the many types defined in Axiom, and produces an output of a particular type. One of the greatest difficulties for the beginner is making sense of the types, and ensuring that your input type is commensurate with the mathematics you are trying to do with it. But this is also one of Axiom’s greatest strengths.

Type Declaration
::
   -> c:PositiiveInteger :=3
       3                       Type:PositiveInteger

In SymPy, to raise something to a power, you must use **, not ^. However, in Axiom, both ^ and ** mean exponentiation.

"""""""""""""
Integration
"""""""""""""
**SymPy**

The integrals module in SymPy implements methdos calculating definite and indefinite integrals of expressions.

Principal method in this module is integrate()
Given below are a few examples
::
     >>> from sympy import *
     >>> import sys
     >>> x=Symbol('X')
     >>> integrate((x**2+2*x+1)/((x+1)**6+1),x)
       atan(X**3 + 3*X**2 + 3*X + 1)/3

Another Example
::
     >>> integrate((sinh(1+sqrt(x+b))+2*sqrt(x+b)) / (sqrt(x+b) * (x + cosh(1+sqrt(x
       + b)))), x)
       2*log(X + cosh((B + X)**(1/2) + 1))
Another Example
::
    >>>integrate(tan(atan(x)/3),x)
    Integral(tan(atan(x)/3),x)

Integrating Error Functions
::
    >>> integrate(exp(-x**2) * erf(x) / (erf(x)**3 - erf(x)**2 - erf(x) + 1),x)
        (erf(x)-1)πlog(erf(x)-1erf(x)+1)-2π8erf(x)-8

**Axiom**

Integration is performed with the integrate command. Sometimes Axiom gets worried that an integrand may have a pole in the region of integration, in which case the addition of the string "noPole" will alleviate its fears:
Here is an example of integration :
::
   -> integrate((x**2+2*x+1)/((x+1)**6+1),x)
       arctan(x^3+3 x^2+3 x+1)/3
                                     Type - Union(Expression Integer,...)

Another Example
::
    -> integrate((sinh(1+sqrt(x+b))+2*sqrt(x+b)) / (sqrt(x+b) * (x + cosh(1+sqrt(x
       + b)))), x)
                              +-----+
                    - 2 cosh(\|x + b  + 1) - 2 x            +-----+
        2 log(---------------------------------------) - 2 \|x + b
                    +-----+              +-----+
              sinh(\|x + b  + 1) - cosh(\|x + b  + 1)
                                     Type: Union(Expression Integer,...)

Another Example

**A strong structure-checking algorithm, Risch Algorithm in Axiom finds hidden algebraic relationships between functions whereas SymPy is unable to integrate it.**

Risch, is an algorithm for the calculus operation of indefinite integration. The algorithm transforms the problem of integration into a problem in algebra. The general case has been solved and implemented in Axiom. This is one of the strongest features of Axiom.
::
    -> integrate(tan(atan(x)/3),x)
        (8 log(3 tan(arctan(x)/3)^2-1)-3 tan(arctan(x)/3)^2+18 x tan(arctan(x)/3))/18
                                    Type: Union(Expression Integer,...)
Integrating Error Functions
::
      integrate(exp(-x**2) * erf(x) / (erf(x)**3 - erf(x)**2 - erf(x) + 1),x)
                     +---+    erf(x) - 1       +---+
        (erf(x) - 1)\|%pi log(----------) - 2 \|%pi
                              erf(x) + 1
        --------------------------------------------
                        8 erf(x) - 8
                                    Type: Union(Expression Integer,...)

However, Sympy can compute a very very difficult integration well Axiom fails:
::
    >>> integrate((x**2+2*x+1+(3*x+1)*sqrt(x+log(x)))/(x*sqrt(x+log(x))*(x+sqrt(x+log(x)))), x)

''''''''''
 MATRIX
''''''''''

**SymPy**

Defining a Matrix

::
  >>> Matrix([[2,-1,0],[-1,2,-1],[0,-1,2]])
      [ 2, -1,  0]
      [-1,  2, -1]
      [ 0, -1,  2]

Inverse of a Matrix
::
   >>> M=Matrix([[2,-1,0],[-1,2,-1],[0,-1,2]])
   >>> M.inv()
   [3/4, 1/2, 1/4]
   [1/2,   1, 1/2]
   [1/4, 1/2, 3/4]

Determinant of a matrix
::
   >>> M.det()
   4

Transpose of a matrix
::
   >>> M.transpose()
   [ 2, -1,  0]
   [-1,  2, -1]
   [ 0, -1,  2]

**Axiom**

Defining a matrix
::
   T:= matrix [[2,-1,0],[-1,2,-1],[0,-1,2]]

        + 2   - 1   0 +
        |             |
        |- 1   2   - 1|
        |             |
        + 0   - 1   2 +
                                  Type: Matrix Integer

The inverse of Matrix

If we look at the inverse of the T matrix we see:
::
        T^-1

        +3  2  1+
        |       |
        |2  2  1|
        |       |
        +1  1  1+
                                   Type: Matrix Fraction Integer


Determinant of a Matrix
::
   1
                                   Type: Fraction Integer


'''''''''''''''''''''''
 Differential Equations
'''''''''''''''''''''''

**Axiom**

An Example
::
   Let y be the unknown function in terms of x.
   -> y := operator 'y
             y
                                 Type: BasicOperator
   -> deq := x**3 * D(y x, x, 3) + x**2 * D(y x, x, 2) - 2 * x * D(y x, x) + 2 * y x = 2 * x**4
         3 ,,,       2 ,,          ,                  4
        x y   (x) + x y  (x) - 2 xy (x) + 2 y(x) = 2 x
                                 Type: Equation Expression Integer

     -> solve(deq, y, x)
                       5       3       2                  3      2       3       3      2
                      x  - 10 x  + 20 x  + 4           2 x  - 3 x  + 1  x  - 1  x  - 3 x  - 1
        [particular = ----------------------, basis = [---------------, ------, -------------]]
                               15 x                           x            x          x
    Type: Union(Record(particular: Expression Integer,basis: List Expression Integer),...)

**SymPy**

An Example
::
   >>>f(x).diff(x, x) + f(x)
   >>>dsolve(f(x).diff(x, x) + f(x), f(x))
   f(x) = C₁⋅sin(x) + C₂⋅cos(x)


'''''''''''''''
Differentiation
'''''''''''''''

**Axiom**

Use the Axiom function 'differentiate' to differentiate an expression.
The differentiate function takes three arguments:
(1) the function to differentiate.
(2) the list of variables for differentiation.
(3) and a list of powers for the variables.

Examples
::
   -> differentiate(cos(x^3),x)
       -3*x**2*sin(x**3)
                                   Type :Expression Integer
   -> differentiate(tan(x),x,1)
       tan(x)^2 + 1
                                   Type :Expression Integer

**SymPy**

Examples
::
   >>> diff(cos(x**3), x)
   -3*X**2*sin(X**3)
   >>> diff(tan(x),x)
   tan(X)**2 + 1
   >>> diff(x**2+x,x)
   2*X + 1

''''''''''''
Conclusion
''''''''''''

SymPy and Axiom are good open source alternatives to Maple/Mathematica. Their goal is to become reasonably fast.
Axiom is a general purpose Computer Algebra system. It is useful for research and development of mathematical algorithms. It defines a strongly typed, mathematically correct type hierarchy. It has a programming language and a built-in compiler.
Efforts are underway to extend this software to :

-  develop a better user interface.
-  make it useful as a teaching tool.
-  develop an algebra server protocol.
-  integrate additional mathematics.
-  rebuild the algebra in a literate programming style.
-  integrate logic programming.
-  develop an Axiom Journal with refereed submissions.

Though Axiom has been acclaimed as a milestone in the field but to date has tended to be used mainly by serious mathematical researchers only.

Another advantage of SymPy is that since it is written in pure Python (and doesn't need anything else), it is perfectly multiplatform, it's small and easy to install and use.

You can choose to use either SymPy or Axiom, depending on what your needs are. For more information you can go to the official sites of SymPy_ and Axiom_.

.. _SymPy: http://sympy.org/
.. _Sage: http://axiom-developer.org/