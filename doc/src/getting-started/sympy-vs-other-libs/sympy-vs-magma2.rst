Differences Between SymPy and Magma
===================================

Brief History
*************

SymPy:
++++++
SymPy is an open source Python library for symbolic mathematics. It only can be used with / in Python. SymPy was created in 2006 and is free to use and to develop. Can be downloaded at `http://sympy.org/en/download.html <http://sympy.org/en/download.html>`_

Magma:
++++++
Magma is a Computer Algebra System (CAS) designed to solve problems in algebra, number theory, geometry and combinatorics. It was officially released in 1993 and created by University of Sydney. Magma is not an open source program. You can order it at `http://magma.maths.usyd.edu.au/magma/ordering/ <http://magma.maths.usyd.edu.au/magma/ordering/>`_

Differences
***********
The first and biggest difference between SymPy and Magma is that SymPy is an open source while Magma is not. It makes SymPy much better than Magma.

The second difference is that SymPy, until now, can be used only with Python language (with .py filetype), while Magma has its own programming language (with .magma filetype).

Magma focuses on Algebra (like Algebraic Number Theory, Linear Algebra, Commutative Algebra, Associative Algebra, etc). Although there are other mathematical areas that can be covered by Magma. Such as Cryptography, Group and Number Theory, Lia Theory, Coding Theory, Optimization, etc. While SymPy focuses on almost all Symbolic Mathematics. Including Algebra (like Magma), basic arithmetic and expression, mathematical function (like sin, cos, tan, zeta, etc), limits, differentiation, calculus, etc. However, there are some weaknesses and advantages in every Mathematical System. And there are differences between areas that are covered by these 2 Mathematical Systems.

The third difference is in the printing way. There are many ways how expression can be printed. I'll give you some example.

**Standard:**

 :math:`x**2 / (x*8)**3`

**Sympy (Standard):**

 :math:`x**2 / (x*8)**3`

You see that "Standard Printing" can be so confusing, if used with complicated formula.

Magma has a little improvements too in this area. It changes '**' symbol with '^' symbol.

**Magma (Standard):**

 :math:`x^2 / (x*8)^3`

But the best way to have an easy-to-read, even for a complicated formula, is with using SymPy's pprint() method.

**SymPy (pprint() method):**

 :math:`x^2`

 :math:`---------`

 :math:`x^3(24)`

Did you see the difference? pprint() method is the simplest, right?

Real World Example
******************

Let's see some examples of how these 2 Mathematical Systems is solving the same problem.

Problem no. 1
+++++++++++++

**Problem:**

 Using Magma and SymPy, solve this formula: (x*2)^3, with using 'x' as a symbol.

**SymPy:**

.. code:: Python

 from sympy import *
 x = Symbol('x')
 pprint((x*2)**3)

**Magma:**

 (Magma cannot solve symbolic mathematics problem)

**Conclusion**

 In this part, SymPy is proved to be a better way to solve a symbolic mathematics problem.

Problem no. 2
+++++++++++++

**Problem:**

 Using Magma and SymPy, create a list with value from 1 to 10. Iterate it, and times each of them.

**SymPy:**

.. code:: Python

 from sympy import *
 x = [1,2,3,4,5,6,7,8,9,10)
 x_count = 0
 for each_x in x:
     x_count += 1
 i = 0
 y = x[i]
 while(i < count-1):
     z = x[i+1]
     zz = y*z
     print(zz)
     y = zz
     i += 1

**Magma:**

.. code:: Python

 x := [1,2,3,4,5,6,7,8,9,10];
 print &*x;

**Conclusion**

 In this problem, Magma is proved to be a better way in operating a list.


Now, you know some examples of both Magma and SymPy. For more examples / tutorials, please visit `http://docs.sympy.org/ <http://docs.sympy.org/>`_ for SymPy, and `http://magma.maths.usyd.edu.au/magma/documentation/ <http://magma.maths.usyd.edu.au/magma/documentation/>`_ for Magma.


Conclusion
**********
We’re sure that everything in this world, has both weaknesses and advantages. And so do Magma and SymPy.

If you want to decide whether using Magma or SymPy, please consider it based on what you really need. I'll give you some advices.

**Use SymPy if:**

-You have experience with Python or another Object Oriented Programming Language (Like Java, C, C#, etc)

-You are a student, and want to learn to develop using Mathematical System for free

-You want to use an Open Source Program


**Use Magma if:**

-You want to focus on developing Algebraic Program. Because Magma has in-depth coverage of all the major branches of algebra, number theory, algebraic geometry and finite incidence geometry

-You are a students from Harvard Mathematics Department

Examine their features, and find what's the best for you.


**I'm Yosi Pramajaya, and this is my article explaining the differences between Magma and SymPy. It's a work in progress for GCI 2011.**




