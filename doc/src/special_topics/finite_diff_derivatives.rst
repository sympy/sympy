===============================================
Finite Difference Approximations to Derivatives
===============================================

Introduction
============

Finite difference approximations to derivatives is quite important in numerical analysis and in
computational physics. In this tutorial we show how to use SymPy to compute  approximations of
varying accuracy. The hope is that these notes could be useful for the practicing researcher who
is developing code in some language and needs to be able to efficiently generate finite difference
formulae for various approximations.

In order to establish notation, we first state that we envision that there exists a continuous function F of a single
variable x, with F having as many derivatives as desired.  We sample x values uniformly at points along the
real line separated by h.  In most cases we want h to be small in some sense.  F(x) may be expanded
about some point `x_{0}` via the usual Taylor series expansion. Let `x = x_{0} + h`. Then the Taylor expansion is

.. math::

   F(x_{0}+h) = F(x_{0}) + \big(\frac{dF}{dx}\big)_{x_{0}} * h +  \frac{1}{2!} \big(\frac{d^{2}F }{dx^{2}}\big)_{x_{0}}* h^2 +
   \frac{1}{3!} \big(\frac{d^{3}F }{dx^{3}}\big)_{x_{0}}* h^3 + ...

In order to simplify the notation, we now define a set of coefficients `c_{n}`, where

.. math::

   c_{n} := \frac{1}{n!} \big(\frac{d^{n}F }{dx^{n}}\big)_{x_{0}}.

So now our series has the form:

.. math::

   F(x_{0}+h) = F(x_{0}) + c_{1} * h +  c_{2}* h^2 + c_{3}* h^3 + ...


In the following we will only use a finite grid of values `x_{i}` with `i` running from `1,...,N` and the corresponding values of our function
F at these grid points denoted by `F_{i}`.  So the problem is how to generate approximate values for the derivatives of F with the constraint that
we use a subset of the finite set of pairs `(x_{i},F_{i})` of size N.

What follows are  manipulations using SymPy to formulate approximations for derivatives of a given order and to assess its accuracy.
First, we use SymPy to derive the approximations by using a rather brute force method frequently covered in introductory treatments. Later we shall make use of other SymPy functions which get the job done with more efficiency.


A Direct Method Using SymPy Matrices
====================================

If we let `x_{0} = x_{i}`, evaluate the series at `x_{i+1}=x_{i}+ h` and truncate all terms above `O(h^1)` we can solve for the single coefficient `c_{1}` and obtain an approximation to the first derivative:

.. math::

	\big(\frac{dF}{dx}\big)_{x_{0}} \approx \frac{F_{i+1} - F_{i}}{h} + O(h)

where the `O(h)` refers to the lowest order term in the series in `h`.  This establishes that the derivative
approximation is of first order accuracy.  Put another way, if we decide that we can only use the two pairs
`(x_{i},F_{i})` and `(x_{i+1},F_{i+1})` we obtain a "first order accurate" derivative.

In addition to `(x_{i},F_{i})` we next use the two points `(x_{i+1},F_{i+1})` and `(x_{i+2},F_{i+2})`.
Then we have two equations:

.. math::
	F_{i+1} = F_{i} + c_{1}* h + \frac{1}{2}*c_{2}*h^2 + \frac{1}{3!}*c_{3}*h^3 + ...
.. math::
	F_{i+2} = F_{i} + c_{1}* (2h) + \frac{1}{2}*c_{2}*(2h)^2 + \frac{1}{3!}*c_{3}*(2h)^3 + ...

If we again want to find the first derivative (`c_{1}`), we can do that by eliminating the term involving `c_{2}` from
the two equations.  We show how to do it using SymPy.


	>>> from __future__ import print_function
	>>> from sympy import *
	>>> x, x0, h = symbols('x, x_0, h')
	>>> Fi, Fip1, Fip2 = symbols('F_{i}, F_{i+1}, F_{i+2}')
	>>> n = 3 # there are the coefficients c_0=Fi, c_1=dF/dx, c_2=d**2F/dx**2
	>>> c = symbols('c:3')
	>>> def P(x, x0, c, n):
	...     return sum( ((1/factorial(i))*c[i] * (x-x0)**i for i in range(n)) )

Vector of right hand sides:

	>>> R = Matrix([[Fi], [Fip1], [Fip2]])

Now we make a matrix consisting of the coefficients
of the c_i in the nth degree polynomial P.

Coefficients of `c_i` evaluated at `x_i`:

	>>> m11 = P(x0 , x0, c, n).diff(c[0])
	>>> m12 = P(x0 , x0, c, n).diff(c[1])
	>>> m13 = P(x0 , x0, c, n).diff(c[2])

Coefficients of `c_i` evaluated at `x_i + h`:

	>>> m21 = P(x0+h, x0, c, n).diff(c[0])
	>>> m22 = P(x0+h, x0, c, n).diff(c[1])
	>>> m23 = P(x0+h, x0, c, n).diff(c[2])

Coefficients of `c_i` evaluated at `x_i + 2*h`:

	>>> m31 = P(x0+2*h, x0, c, n).diff(c[0])
	>>> m32 = P(x0+2*h, x0, c, n).diff(c[1])
	>>> m33 = P(x0+2*h, x0, c, n).diff(c[2])

Matrix of the coeffcients is 3x3 in this case:

	>>> M = Matrix([[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]])

Matrix form of the three equations for the `c_i` is M*X = R:

The solution is obtained by directly inverting the 3x3 matrix M:

	>>> X =  M.inv() * R

Note that all three coefficients make up the solution. The desired first derivative is coefficient `c_1` which is X[1].

	>>> print(together(X[1]))
	(4*F_{i+1} - F_{i+2} - 3*F_{i})/(2*h)

It is instructive to compute another three-point approximation to the first derivative,  except centering the approximation
at `x_i` and thus using points at `x_{i-1}`,  `x_{i}`,  and `x_{i+1}`. So here is how this can be done using the 'brute force' method:


	>>> from __future__ import print_function
	>>> from sympy import *
	>>> x, x0, h = symbols('x, x_i, h')
	>>> Fi, Fim1, Fip1 = symbols('F_{i}, F_{i-1}, F_{i+1}')
	>>> n = 3 # there are the coefficients c_0=Fi,  c_1=dF/h,  c_2=d**2F/h**2
	>>> c = symbols('c:3')
	>>> # define a polynomial of degree n
	>>> def P(x, x0, c, n):
	...    return sum( ((1/factorial(i))*c[i] * (x-x0)**i for i in range(n)) )
	>>> # now we make a matrix consisting of the coefficients
	>>> # of the c_i in the nth degree polynomial P
	>>> # coefficients of c_i evaluated at x_i
	>>> m11 = P(x0 , x0, c, n).diff(c[0])
	>>> m12 = P(x0 , x0, c, n).diff(c[1])
	>>> m13 = P(x0 , x0, c, n).diff(c[2])
	>>> # coefficients of c_i evaluated at x_i - h
	>>> m21 = P(x0-h, x0, c, n).diff(c[0])
	>>> m22 = P(x0-h, x0, c, n).diff(c[1])
	>>> m23 = P(x0-h, x0, c, n).diff(c[2])
	>>> # coefficients of c_i evaluated at x_i + h
	>>> m31 = P(x0+h, x0, c, n).diff(c[0])
	>>> m32 = P(x0+h, x0, c, n).diff(c[1])
	>>> m33 = P(x0+h, x0, c, n).diff(c[2])
	>>> # matrix of the coeffcients is 3x3 in this case
	>>> M = Matrix([[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]])

Now that we have the matrix of coefficients we next form the right-hand-side and solve by inverting `M`:

	>>> # matrix of the function values...actually a vector of right hand sides
	>>> R = Matrix([[Fi], [Fim1], [Fip1]])
	>>> # matrix form of the three equations for the c_i is M*X = R
	>>> # solution directly inverting the 3x3 matrix M:
	>>> X =  M.inv() * R
	>>> # note that all three coefficients make up the solution
	>>> # the first derivative is coefficient c_1 which is X[1].
	>>> print("The second-order accurate approximation for the first derivative is: ")
	The second-order accurate approximation for the first derivative is:
	>>> print(together(X[1]))
	(F_{i+1} - F_{i-1})/(2*h)

These two examples serve to show how one can directly find second order accurate first derivatives using SymPy.
The first example uses values of `x` and `F` at all three points `x_i`, `x_{i+1}`, and `x_{i+2}` whereas the
second example only uses values of `x` at the two points `x_{i-1}` and `x_{i+1}` and thus is a bit more efficient.

From these two simple examples a general rule is that if one wants a first derivative to be accurate to `O(h^{n})`
then one needs n+1 function values in the approximating polynomial (here provided via the function `P(x,x0,c,n)`).


Now let's assess the question of the accuracy of the centered difference result to see how we determine that it is
really second order.  To do this we take the result for `dF/dx` and substitute in the polynomial expansion for a higher
order polynomial and see what we get. To this end,  we make a set of eight coefficients d and use them to perform the
check:


    >>> d = symbols('c:8')
    >>> dfdxcheck = (P(x0+h, x0, d, 8) - P(x0-h, x0, d, 8))/(2*h)
    >>> print(simplify(dfdxcheck)) # so the appropriate cancellation of terms involving `h` happens
    c1 + c3*h**2/6 + c5*h**4/120 + c7*h**6/5040

Thus we see that indeed the derivative is `c_1` with the next term in the series of order `h^2`.

However,  it can quickly become rather tedious to generalize the direct method as presented above when attempting
to generate a derivative approximation to high order,  such as 6 or 8 although the method certainly works and using
the present method is certainly less tedious than performing the calculations by hand.

As we have seen in the discussion above,  the simple centered approximation for the first derivative only uses two
point values of the `(x_{i},F_{i})` pairs.  This works fine until one encounters the last point in the domain,  say at
`i=N`. Since our centered derivative approximation would use data at the point `(x_{N+1},F_{N+1})` we see that the
derivative formula will not work. So,  what to do?  Well,  a simple way to handle this is to devise a different formula
for this last point which uses points for which we do have values. This is the so-called backward difference formula.
To obtain it,  we can use the same direct approach,  except now us the three points `(x_{N},F_{N})`,  `(x_{N-1},F_{N-1})`,
and `(x_{N-2},F_{N-2})` and center the approximation at `(x_{N},F_{N})`. Here is how it can be done using SymPy:


	>>> from __future__ import print_function
	>>> from sympy import *
	>>> x, xN, h = symbols('x, x_N, h')
	>>> FN, FNm1, FNm2 = symbols('F_{N}, F_{N-1}, F_{N-2}')
	>>> n = 8 # there are the coefficients c_0=Fi,  c_1=dF/h,  c_2=d**2F/h**2
	>>> c = symbols('c:8')
	>>> # define a polynomial of degree d
	>>> def P(x, x0, c, n):
	...     return sum( ((1/factorial(i))*c[i] * (x-x0)**i for i in range(n)) )

Now we make a matrix consisting of the coefficients of the `c_i` in the dth
degree polynomial P coefficients of `c_i` evaluated at `x_i, x_{i-1},` and `x_{i+1}`:

    >>> m11 = P(xN , xN, c, n).diff(c[0])
    >>> m12 = P(xN, xN, c, n).diff(c[1])
    >>> m13 = P(xN , xN, c, n).diff(c[2])
    >>> # coefficients of c_i evaluated at x_i - h
    >>> m21 = P(xN-h, xN, c, n).diff(c[0])
    >>> m22 = P(xN-h, xN, c, n).diff(c[1])
    >>> m23 = P(xN-h, xN, c, n).diff(c[2])
    >>> # coefficients of c_i evaluated at x_i + h
    >>> m31 = P(xN-2*h, xN, c, n).diff(c[0])
    >>> m32 = P(xN-2*h, xN, c, n).diff(c[1])
    >>> m33 = P(xN-2*h, xN, c, n).diff(c[2])

Next we construct the `3 \times 3` matrix of the coeffcients:

    >>> M = Matrix([[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]])
    >>> # matrix of the function values...actually a vector of right hand sides
    >>> R = Matrix([[FN], [FNm1], [FNm2]])

Then we invert `M` and write the solution to the `3 \times 3` system.

The matrix form of the three equations for the c_i is `M*C = R`. The solution is obtained by
directly inverting `M`:

    >>> X =  M.inv() * R

The first derivative is coefficient `c_1` which is `X[1]`. Thus the second order accurate
approximation for the first derivative is:

    >>> print("The first derivative centered at the last point on the right is:")
    The first derivative centered at the last point on the right is:
    >>> print(together(X[1]))
    (-4*F_{N-1} + F_{N-2} + 3*F_{N})/(2*h)

Of course,  we can devise a similar formula for the value of the derivative at the left end
of the set of points at `(x_{1},F_{1})` in terms of values at `(x_{2},F_{2})` and `(x_{3},F_{3})`.

Also,  we note that output of formats appropriate to Fortran,  C,  etc. may be done in the examples
given above.

Next we show how to perform these and many other discritizations of derivatives,  but using a
much more efficient approach originally due to Bengt Fornberg and now incorported into SymPy.



