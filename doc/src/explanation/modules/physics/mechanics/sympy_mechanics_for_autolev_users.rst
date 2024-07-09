.. _sympy_mechanics_for_autolev_users:

=================================
SymPy Mechanics for Autolev Users
=================================

Introduction
------------

Autolev (now superseded by MotionGenesis) is a domain specific programming
language which is used for symbolic multibody dynamics. The SymPy mechanics
module now has enough power and functionality to be a fully featured symbolic
dynamics module. The PyDy package extends the SymPy output to the numerical
domain for simulation, analyses and visualization. Autolev and SymPy Mechanics have
a lot in common but there are also many differences between them.
This page shall expand upon their differences. It is meant to be a go-to
reference for Autolev users who want to transition to SymPy Mechanics.

It would be nice to have a basic understanding of SymPy and SymPy Mechanics before
going over this page.
If you are completely new to Python, you can check out the official
`Python Tutorial <https://docs.python.org/3/tutorial/>`_.
Check out the :ref:`SymPy Documentation <documentation>`, especially
the tutorial to get a feel for SymPy.
For an introduction to Multibody dynamics in Python, `this <https://www.youtube.com/watch?v=mdo2NYtA-xY&t=6950s>`_
lecture is very helpful.

You might also find the :ref:`Autolev Parser <autolev_parser>` which is
a part of SymPy to be helpful.

Some Key Differences
------------------------

+-----------------------------------+-----------------------------------+
|          **Autolev**              |         **SymPy Mechanics**       |
+===================================+===================================+
||                                  ||                                  |
| Autolev is a domain specific      | SymPy is a library written in the |
| programming language designed to  | general purpose language Python.  |
| perform multibody dynamics. Since | Although Autolev's code is more   |
| it is a language of its own, it   | compact, SymPy (by virtue of being|
| has a very rigid language         | an add on to Python) is more      |
| specification. It predefines,     | flexible. The users have more     |
| assumes and computes              | control over what they can do. For|
| many things based on the          | example, one can create a class in|
| input code. Its code is a lot     | their code for let's say a type of|
| cleaner and concise as a result of| rigibodies with common            |
| this.                             | properties.                       |
|                                   | The wide array of scientific      |
|                                   | Python libraries available is also|
|                                   | a big plus.                       |
+-----------------------------------+-----------------------------------+
||                                  ||                                  |
| Autolev generates Matlab, C, or   | SymPy generates numerical Python, |
| Fortran code from a small set of  | C or Octave/Matlab code from a    |
| symbolic mathematics.             | large set of symbolic mathematics |
|                                   | created with SymPy. It also builds|
|                                   | on the popular scientific Python  |
|                                   | stack such as NumPy, SciPy,       |
|                                   | IPython, matplotlib, Cython and   |
|                                   | Theano.                           |
+-----------------------------------+-----------------------------------+
||                                  ||                                  |
| Autolev uses 1 (one) based        | Python uses 0 (zero) based        |
| indexing. The initial element of  | indexing. The initial element of  |
| a sequence is found using a[1].   | a sequence is found using a[0].   |
+-----------------------------------+-----------------------------------+
||                                  ||                                  |
| Autolev is case insensitive.      | SymPy code being Python code is   |
|                                   | case sensitive.                   |
+-----------------------------------+-----------------------------------+
||                                  ||                                  |
| One can define their own commands | SymPy code is Python code, so one |
| in Autolev by making .R and .A    | can define functions in their     |
| files which can be used in their  | code. This is a lot more          |
| programs.                         | convenient.                       |
+-----------------------------------+-----------------------------------+
||                                  ||                                  |
| Autolev is proprietary.           | SymPy is open source.             |
+-----------------------------------+-----------------------------------+

Rough Autolev-SymPy Equivalents
----------------------------------

The tables below give rough equivalents for some common Autolev
expressions. **These are not exact equivalents**, but rather should be
taken as hints to get you going in the right direction. For more detail
read the built-in documentation on :ref:`SymPy vectors <physics_vector>`,
:ref:`SymPy mechanics <physics_mechanics>` and
`PyDy <https://www.pydy.org/documentation.html>`_ .

In the tables below, it is assumed that you have executed the following
commands in Python:
::

    import sympy.physics.mechanics as me
    import sympy as sm

Mathematical Equivalents
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+-----------------------+-----------------------+-----------------------+
| **Autolev**           | **SymPy**             | **Notes**             |
+=======================+=======================+=======================+
||                      ||                      ||                      |
| ``Constants A, B``    | ``a, b =              | Note that the names   |
|                       | sm.symbols(‘a         | of the symbols can be |
|                       |  b’, real=True)``     | different from the    |
|                       |                       | names of the          |
|                       |                       | variables they are    |
|                       |                       | assigned to. We can   |
|                       |                       | define ``a, b =       |
|                       |                       | symbols(‘b a’)`` but  |
|                       |                       | its good practice to  |
|                       |                       | follow the            |
|                       |                       | convention.           |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| ``Constants C+``      | ``c = sm.symbols(‘c’, | Refer to SymPy        |
|                       | real=True,            | :ref:`assumptions     |
|                       | nonnegative=True)``   | <assumptions_module>` |
|                       |                       | for more information. |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| ``Constants D-``      | ``d = sm.symbols(‘d’, |                       |
|                       | real=True,            |                       |
|                       | nonpositive=True)``   |                       |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| ``Constants K{4}``    | ``k1, k2, k3, k4 =    |                       |
|                       | sm.symbols('k1 k2 k3  |                       |
|                       | k4', real=True)``     |                       |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| ``Constants a{2:4}``  | ``a2, a3, a4 =        |                       |
|                       | sm.symbols('a2 a3 a4',|                       |
|                       | real=True)``          |                       |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| ``Constants           | ``b11, b12, b21, b22 =|                       |
| b{1:2, 1:2}``         | sm.symbols('b11 b12   |                       |
|                       | b21 b22', real=True)``|                       |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| ``Specified Phi``     | ``phi =               |                       |
|                       | me.dynamicsymbols(‘phi|                       |
|                       | ')``                  |                       |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| ``Variables q, s``    | ``q, s =              |                       |
|                       | me.dynamicsymbols(q,  |                       |
|                       | s)``                  |                       |
+-----------------------+-----------------------+-----------------------+
| ``Variables x’’``     | ``x =                 |                       |
|                       | me.dynamicsymbols(‘x’ |                       |
|                       | )``                   |                       |
|                       |                       |                       |
|                       | ``xd =                |                       |
|                       | me.dynamicsymbols(‘x’ |                       |
|                       | , 1)``                |                       |
|                       |                       |                       |
|                       | ``xd2 =               |                       |
|                       | me.dynamicsymbols(‘x’ |                       |
|                       | , 2)``                |                       |
+-----------------------+-----------------------+-----------------------+
| ``Variables y{2}’``   | ``y1 =                |                       |
|                       | me.dynamicsymbols(‘y1’|                       |
|                       | )``                   |                       |
|                       |                       |                       |
|                       | ``y2 =                |                       |
|                       | me.dynamicsymbols(‘y2’|                       |
|                       | )``                   |                       |
|                       |                       |                       |
|                       | ``y1d =               |                       |
|                       | me.dynamicsymbols(‘y1’|                       |
|                       | , 1)``                |                       |
|                       |                       |                       |
|                       | ``y2d =               |                       |
|                       | me.dynamicsymbols(‘y2'|                       |
|                       | , 1)``                |                       |
+-----------------------+-----------------------+-----------------------+
| ``MotionVariables     | ``u1 =                | SymPy doesn’t         |
| u{2}``                | me.dynamicsymbols(‘u1’| differentiate between |
|                       | )``                   | variables,            |
|                       |                       | motionvariables and   |
|                       | ``u2 =                | specifieds during     |
|                       | me.dynamicsymbols('u2'| declaration. Instead, |
|                       | )``                   | it takes different    |
|                       |                       | lists of these as     |
|                       |                       | parameters in objects |
|                       |                       | like the KanesMethod. |
+-----------------------+-----------------------+-----------------------+
| ``Imaginary j``       | ``j = sm.I``          | I is a sympy object   |
|                       |                       | which stands for the  |
|                       |                       | imaginary unit. One   |
|                       |                       | can define complex    |
|                       |                       | numbers using it.     |
|                       |                       |                       |
|                       |                       | ``z = x + I*y``       |
|                       |                       |                       |
|                       |                       | where x, y and z are  |
|                       |                       | symbols.              |
+-----------------------+-----------------------+-----------------------+
| ``Tina = 2*pi``       | ``tina = 2*sm.pi``    | Using ``.evalf()``    |
|                       |                       | will result in the    |
|                       | ``tina =              | numeric value.        |
|                       | tina.evalf()``        |                       |
|                       |                       |                       |
| ``s = u*t + a*t^2/2`` | ``t =                 |                       |
|                       | me.dynamicsymbols._t``|                       |
|                       |                       |                       |
|                       | ``s = u*t + a*t**2/2``|                       |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| ``abs(x)^3 + sin(x)^2 | ``sm.abs(x)**3        |                       |
| + acos(x)``           | + sm.sin(x)**2 +      |                       |
|                       | sm.acos(x)``          |                       |
+-----------------------+-----------------------+-----------------------+
| ``E = (x+2*y)^2 +     | ``E = (x+2*y)**2 +    | For more information  |
| 3*(7+x)*(x+y)``       | 3*(7+x)*(x+y)``       | refer to              |
|                       |                       | :ref:`simplification. |
| ``Expand(E)``         | ``sm.expand(E)``      | <tutorial-simplify>`  |
|                       |                       |                       |
| ``Factor(E, x)``      | ``sm.horner(E,        |                       |
|                       | wrt=x)``              |                       |
|                       |                       |                       |
| ``Coef(y, x)``        | ``y.coeff(x)``        | These SymPy functions |
|                       |                       | do not work in place. |
| ``Replace(y,          | ``y.subs({sm.sin(x):  | They just return      |
| sin(x)=3)``           | 3})``                 | expressions. If you   |
|                       |                       | want to overwrite the |
| ``Exclude(E,x)``      | ``e.collect(x).coeff( | original expression   |
|                       | x, 0)``               | you would have to do  |
|                       |                       | something like:       |
| ``Include(E,x)``      | ``e.collect(x).coeff( |                       |
|                       | x, 1)``               | ``y =                 |
|                       |                       | y.subs({sm.sin(x):    |
| ``Arrange(E,2,y)``    | ``e.collect(y)``      | 3})``                 |
+-----------------------+-----------------------+-----------------------+
| ``Dy = D(E, y)``      | ``E.diff(y)``         | For more information  |
|                       |                       | refer to              |
| ``Dt = Dt(E)``        | ``E.diff(             | :ref:`calculus.       |
|                       | me.dynamicsymbols._t  | <calculus>`           |
|                       | )``                   |                       |
|                       |                       |                       |
|                       | Works if the          |                       |
|                       | expression is made up |                       |
|                       | of dynamicsymbols.    |                       |
| ``Dt2 = Dt(V, A)``    |                       |                       |
| where V is a vector   | ``dt2 = v.dt(A)``     |                       |
| and A is a frame      |                       |                       |
|                       |                       |                       |
| ``Dy2 = D(V, y, A)``  | ``dy2 = v.diff(y, A)``|                       |
|                       |                       |                       |
+-----------------------+-----------------------+-----------------------+
| ``E = COS(X*Y)``      | ``e = sm.cos(x*y)``   | For more information  |
|                       |                       | refer to :ref:`series.|
| ``TY = Taylor(E,      | ``b = e.series(x, 0,  | <series_expansions>`  |
| 0:2, x=0, y=0)``      | 2).removeO().series(y,|                       |
|                       | 0, 2).removeO()``     |                       |
|                       |                       |                       |
+-----------------------+-----------------------+-----------------------+
| ``F = Evaluate(E,     | ``E.subs([(x, a), (y, |                       |
| x=a, y=2)``           | 2)])``                |                       |
|                       |                       |                       |
|                       | To get floating point |                       |
|                       | numbers from numerical|                       |
|                       | expressions use       |                       |
|                       | ``.evalf()``          |                       |
|                       |                       |                       |
|                       | ``E.evalf((a +        |                       |
|                       | sm.pi).subs({a: 3}))``|                       |
+-----------------------+-----------------------+-----------------------+
| ``P = Polynomial([a,  | ``p =                 | For more information  |
| b, c], x)``           | sm.Poly(sm.Matrix([a, | refer to              |
|                       | b, c]).reshape(1, 3), | :ref:`polys.          |
|                       | x)``                  | <polys-reference>`    |
+-----------------------+-----------------------+-----------------------+
| ``Roots(Polynomial(   | ``sm.solve(           | For more information  |
| a*x^2 + b*x + c, x,   | sm.Poly(a*x**2 +      | refer to              |
| 2)``                  | b*x + c))``           | :ref:`solvers-docs`.  |
|                       |                       |                       |
| ``Roots([1;2;3])``    | ``sm.solve(sm.Poly(   | For numerical         |
|                       | sm.Matrix([1,2,3]).   | computation related   |
|                       | reshape(3, 1), x),    | to polynomials and    |
|                       | x)``                  | roots refer to        |
|                       |                       | `mpmath/calculus. <ht |
|                       |                       | tps://web.archive.org |
|                       |                       | /web/20180731093609/h |
|                       |                       | ttp://docs.sympy.org/ |
|                       |                       | 0.7.6/modules/mpmath/ |
|                       |                       | calculus/polynomials. |
|                       |                       | html>`_               |
+-----------------------+-----------------------+-----------------------+
| ``Solve(A, x1, x2)``  | ``sm.linsolve(A,      | For more information  |
|                       | (x1, x2))``           | refer to              |
|                       |                       | :ref:`                |
| where A is an         | where A is an         | solvers/solveset.     |
| augmented matrix that | augmented matrix      | <solveset>`           |
| represents the linear |                       |                       |
| equations and x1, x2  |                       |                       |
| are the variables to  |                       | For non linear solvers|
| solve for.            |                       | refer to              |
|                       |                       | ``nonlinsolve`` and   |
|                       |                       | ``nsolve`` in         |
|                       |                       | :ref:`solvers.        |
|                       |                       | <solvers-docs>`       |
+-----------------------+-----------------------+-----------------------+
| ``RowMatrix = [1, 2,  | ``row_matrix =        | For more information  |
| 3, 4]``               | sm.Matrix([[1],[2],   | refer to              |
|                       | [3],[4]])``           | :ref:`matrices.       |
|                       |                       | <matrices>`           |
| ``ColMatrix = [1; 2;  | ``col_matrix =        |                       |
| 3; 4]``               | sm.Matrix([1, 2, 3,   |                       |
|                       | 4])``                 |                       |
|                       |                       |                       |
| ``MO = [a, b; c, 0]`` | ``MO = sm.Matrix([[a, |                       |
|                       | b], [c, 0]])``        |                       |
|                       |                       |                       |
| ``MO[2, 2] := d``     | ``MO[1, 1] = d``      |                       |
|                       |                       |                       |
| ``A + B*C``           | ``A + B*C``           |                       |
|                       |                       |                       |
| ``Cols(A)``           | ``A.cols``            |                       |
|                       |                       |                       |
| ``Cols(A, 1)``        | ``A.col(0)``          |                       |
|                       |                       |                       |
| ``Rows(A)``           | ``A.rows``            |                       |
|                       |                       |                       |
| ``Rows(A, 1)``        | ``A.row(0)``          |                       |
|                       |                       |                       |
| ``Det(A)``            | ``M.det()``           |                       |
|                       |                       |                       |
| ``Element(A, 2, 3)``  | ``M[2, 3]``           |                       |
|                       |                       |                       |
| ``Inv(A)``            | ``M**-1``             |                       |
|                       |                       |                       |
| ``Trace(A)``          | ``sm.trace(A)``       |                       |
|                       |                       |                       |
| ``Transpose(A)``      | ``A.T``               |                       |
|                       |                       |                       |
| ``Diagmat(4, 1)``     | ``sm.diag(1,1,1,1)``  |                       |
|                       |                       |                       |
| ``Eig(A)``            | ``A.eigenvals()``     |                       |
|                       |                       |                       |
| ``Eig(A, EigVal,      | ``eigval =            |                       |
| EigVec)``             | A.eigenvals()``       |                       |
|                       |                       |                       |
|                       | ``eigvec =            |                       |
|                       | A.eigenvects()``      |                       |
+-----------------------+-----------------------+-----------------------+


Physical Equivalents
~~~~~~~~~~~~~~~~~~~~~~~~

+-----------------------+-----------------------+-----------------------+
| **Autolev**           | **SymPy**             | **Notes**             |
+=======================+=======================+=======================+
| ``Bodies A``          | ``m =sm.symbols(‘m’)``| The 4th and 5th       |
|                       |                       | arguments are for the |
| Declares A, its       | ``Ao =                | mass and inertia.     |
| masscenter Ao, and    | sm.symbols(‘Ao’)``    | These are specified   |
| orthonormal vectors   |                       | after the declaration |
| A1>, A2> and A3>      | ``Af =                | in Autolev.           |
| fixed in A.           | me.ReferenceFrame(‘Af’|                       |
|                       | )``                   |                       |
|                       |                       |                       |
|                       | ``I =                 | One can pass a dummy  |
|                       | me.outer(Af.x,Af.x)`` | for the parameters    |
|                       |                       | and use setters       |
|                       | ``P = me.Point(‘P’)`` | ``A.mass = \_`` and   |
|                       |                       | ``A.inertia = \_`` to |
|                       | ``A =me.RigidBody(‘A’,| set them later.       |
|                       | Ao, Af, m, (I, P))``  |                       |
|                       |                       |                       |
|                       | Af.x, Af.y and Af.z   | For more information  |
|                       | are equivalent to     | refer to              |
|                       | A1>, A2> and A3>.     | :ref:`mechanics/masses|
|                       |                       | .<masses>`            |
+-----------------------+-----------------------+-----------------------+
| ``Frames A``          | ``A =                 | For more information  |
|                       | me.ReferenceFrame(‘A’ | refer to              |
| ``V1> =               | )``                   | :ref:`physics/vectors.|
| X1*A1> + X2*A2>``     |                       | <matrices>`           |
|                       | ``v1 =                |                       |
|                       | x1*A.x + x2*A.y``     |                       |
+-----------------------+-----------------------+-----------------------+
| ``Newtonian N``       | ``N =                 | SymPy doesn’t specify |
|                       | me.ReferenceFrame(‘N’ | that a frame is       |
|                       | )``                   | inertial during       |
|                       |                       | declaration. Many     |
|                       |                       | functions such as     |
|                       |                       | ``set_ang_vel()`` take|
|                       |                       | the inertial          |
|                       |                       | reference frame as a  |
|                       |                       | parameter.            |
+-----------------------+-----------------------+-----------------------+
| ``Particles C``       | ``m =                 | The 2nd and 3rd       |
|                       | sm.symbols(‘m’)``     | arguments are for the |
|                       |                       | point and mass. In    |
|                       | ``Po =                | Autolev, these are    |
|                       | me.Point(‘Po’)``      | specified after the   |
|                       |                       | declaration.          |
|                       | ``C = me.Particle(‘C’,|                       |
|                       | Po, m)``              | One can pass a dummy  |
|                       |                       | and use setters       |
|                       |                       | (``A.point = \_`` and |
|                       |                       | ``A.mass = \_``) to   |
|                       |                       | set them later.       |
+-----------------------+-----------------------+-----------------------+
| ``Points P, Q``       | ``P = me.Point(‘P’)`` |                       |
|                       |                       |                       |
|                       | ``Q = me.Point(‘Q’)`` |                       |
+-----------------------+-----------------------+-----------------------+
| ``Mass B=mB``         | ``mB = symbols(‘mB’)``|                       |
|                       |                       |                       |
|                       | ``B.mass = mB``       |                       |
+-----------------------+-----------------------+-----------------------+
| ``Inertia B, I1, I2,  | ``I = me.inertia(Bf,  | For more information  |
| I3, I12, I23, I31``   | i1, i2, i3, i12, i23, | refer to the          |
|                       | i31)``                | :ref:`mechanics api.  |
|                       |                       | <part_bod>`           |
|                       | ``B.inertia = (I, P)``|                       |
|                       | where B is a          |                       |
|                       | rigidbody, Bf is the  |                       |
|                       | related frame and P is|                       |
|                       | the center of mass of |                       |
|                       | B.                    |                       |
|                       |                       |                       |
|                       | Inertia dyadics can   |                       |
|                       | also be formed using  |                       |
|                       | vector outer products.|                       |
|                       |                       |                       |
|                       | ``I =                 |                       |
|                       | me.outer(N.x, N.x)``  |                       |
+-----------------------+-----------------------+-----------------------+
| ``vec> = P_O_Q>/L``   | ``vec  =              | For more information  |
|                       | (Qo.pos_from(O))/L``  | refer to              |
| ``vec> =              |                       | :ref:`physics/vectors.|
| u1*N1> + u2*N2>``     | ``vec =               | <physics_vector>`     |
|                       | u1*N.x + u2*N.y``     |                       |
| ``Cross(a>, b>)``     |                       |                       |
|                       | ``cross(a, b)``       |                       |
| ``Dot(a>, b>)``       |                       |                       |
|                       | ``dot(a, b)``         |                       |
| ``Mag(v>)``           |                       |                       |
|                       | ``v.magnitude()``     |                       |
| ``Unitvec(v>)``       |                       |                       |
|                       | ``v.normalize()``     |                       |
|                       |                       |                       |
| ``DYAD>> = 3*A1>*A1> +| ``dyad =              |                       |
| A2>*A2> + 2*A3>*A3>`` | 3*me.outer(a.x        |                       |
|                       | ,a.x) + me.outer(a.y, |                       |
|                       | a.y) + 2*me.outer(a.z |                       |
|                       | ,a.z)``               |                       |
+-----------------------+-----------------------+-----------------------+
| ``P_O_Q> = LA*A1>``   | ``Q.point =           | For more information  |
|                       | O.locatenew(‘Qo’,     | refer to the          |
|                       | LA*A.x)``             | :ref:`kinematics api. |
|                       |                       | <kinematics>`         |
| ``P_P_Q> = LA*A1>``   | where A is a          |                       |
|                       | reference frame.      |                       |
|                       |                       |                       |
|                       | ``Q.point =           |                       |
|                       | P.point.locatenew(‘Qo | All these vector and  |
|                       | ’,                    | kinematic functions   |
|                       | LA*A.x)``             | are to be used on     |
|                       |                       | ``Point`` objects and |
|                       |                       | not ``Particle``      |
|                       |                       | objects so ``.point`` |
|                       |                       | must be used for      |
|                       |                       | particles.            |
+-----------------------+-----------------------+-----------------------+
| ``V_O_N> = u3*N.1> +  | ``O.set_vel(N, u1*N.x | The getter would be   |
| u4*N.2>``             | + u2*N.y)``           | ``O.vel(N)``.         |
|                       |                       |                       |
| ``Partials(V_O_N>,    | ``O.partial_velocity(N|                       |
| u3)``                 | , u3)``               |                       |
+-----------------------+-----------------------+-----------------------+
| ``A_O_N> = 0>``       | ``O.set_acc(N, 0)``   | The getter would be   |
|                       |                       | ``O.acc(N)``.         |
| Acceleration of point |                       |                       |
| O in reference frame  |                       |                       |
| N.                    |                       |                       |
+-----------------------+-----------------------+-----------------------+
| ``W_B_N> = qB’*B3>``  | ``B.set_ang_vel(N,    | The getter would be   |
|                       | qBd*Bf.z)``           | ``B.ang_vel_in(N)``.  |
| Angular velocity of   |                       |                       |
| body B in reference   | where Bf is the frame |                       |
| frame F.              | associated with the   |                       |
|                       | body B.               |                       |
+-----------------------+-----------------------+-----------------------+
| ``ALF_B_N> =Dt(W_B_N>,| ``B.set_ang_acc(N,    | The getter would be   |
| N)``                  | diff(B.ang_vel_in(N)  | ``B.ang_acc_in(N)``.  |
|                       | )``                   |                       |
| Angular acceleration  |                       |                       |
| of body B in          |                       |                       |
| reference frame N.    |                       |                       |
+-----------------------+-----------------------+-----------------------+
| ``Force_O> = F1*N1> + | In SymPy one should   |                       |
| F2*N2>``              | have a list which     |                       |
|                       | contains all the      |                       |
| ``Torque_A> =         | forces and torques.   |                       |
| -c*qA’*A3>``          |                       |                       |
|                       | ``fL.append((O, f1*N.x|                       |
|                       | + f2*N.y))``          |                       |
|                       |                       |                       |
|                       | where fL is the force |                       |
|                       | list.                 |                       |
|                       |                       |                       |
|                       | ``fl.append((A,       |                       |
|                       | -c*qAd*A.z))``        |                       |
+-----------------------+-----------------------+-----------------------+
| ``A_B = M``           | ``B.orient(A, 'DCM',  |                       |
| where M is a matrix   | M)`` where M is a     |                       |
| and A, B are frames.  | SymPy Matrix.         |                       |
|                       |                       |                       |
| ``D = A_B*2 + 1``     | ``D = A.dcm(B)*2 + 1``|                       |
+-----------------------+-----------------------+-----------------------+
| ``CM(B)``             | ``B.masscenter``      |                       |
+-----------------------+-----------------------+-----------------------+
| ``Mass(A,B,C)``       | ``A.mass + B.mass +   |                       |
|                       | C.mass``              |                       |
+-----------------------+-----------------------+-----------------------+
| ``V1pt(A,B,P,Q)``     | ``Q.v1pt_theory(P, A, | P and Q are assumed to|
|                       | B)``                  | be ``Point`` objects  |
|                       |                       | here. Remember to use |
|                       |                       | ``.point`` for        |
|                       |                       | particles.            |
+-----------------------+-----------------------+-----------------------+
| ``V2pts(A,B,P,Q)``    | ``Q.v2pt_theory(P, A, |                       |
|                       | B)``                  |                       |
+-----------------------+-----------------------+-----------------------+
| ``A1pt(A,B,P,Q)``     | ``Q.a1pt_theory(P, A, |                       |
|                       | B)``                  |                       |
+-----------------------+-----------------------+-----------------------+
| ``A2pts(A,B,P,Q)``    | ``Q.a2pt_theory(P, A, |                       |
|                       | B)``                  |                       |
+-----------------------+-----------------------+-----------------------+
| ``Angvel(A,B)``       | ``B.ang_vel_in(A)``   |                       |
+-----------------------+-----------------------+-----------------------+
| ``Simprot(A, B, 1,    | ``B.orient(A, ‘Axis’, |                       |
| qA)``                 | qA, A.x)``            |                       |
+-----------------------+-----------------------+-----------------------+
| ``Gravity(G*N1>)``    | ``fL.extend(gravity(  | In SymPy we must use a|
|                       | g*N.x, P1, P2, ...))``| forceList (here fL)   |
|                       |                       | which contains tuples |
|                       |                       | of the form ``(point, |
|                       |                       | force_vector)``. This |
|                       |                       | is passed to the      |
|                       |                       | ``kanes_equations()`` |
|                       |                       | method of the         |
|                       |                       | KanesMethod object.   |
+-----------------------+-----------------------+-----------------------+
| ``CM(O,P1,R)``        | ``me.functions.       |                       |
|                       | center_of_mass(o, p1, |                       |
|                       | r)``                  |                       |
+-----------------------+-----------------------+-----------------------+
| ``Force(P/Q, v>)``    | ``fL.append((P, -1*v),|                       |
|                       | (Q, v))``             |                       |
+-----------------------+-----------------------+-----------------------+
| ``Torque(A/B, v>)``   | ``fL.append((A, -1*v),|                       |
|                       | (B, v))``             |                       |
+-----------------------+-----------------------+-----------------------+
| ``Kindiffs(A, B ...)``| ``KM.kindiffdict()``  |                       |
+-----------------------+-----------------------+-----------------------+
| ``Momentum(option)``  | ``linear_momentum(N,  |                       |
|                       | B1, B2 ...)``         |                       |
|                       |                       |                       |
|                       | reference frame       |                       |
|                       | followed by one or    |                       |
|                       | more bodies           |                       |
|                       |                       |                       |
|                       | ``angular_momentum(O, |                       |
|                       | N, B1, B2 ...)``      |                       |
|                       |                       |                       |
|                       | point, reference      |                       |
|                       | frame followed by one |                       |
|                       | or more bodies        |                       |
+-----------------------+-----------------------+-----------------------+
| ``KE()``              | ``kinetic_energy(N,   |                       |
|                       | B1, B2 ...)``         |                       |
|                       |                       |                       |
|                       | reference frame       |                       |
|                       | followed by one or    |                       |
|                       | more bodies           |                       |
+-----------------------+-----------------------+-----------------------+
| ``Constrain(...)``    | ``velocity_constraints| For more details      |
|                       | = [...]``             | refer to              |
|                       |                       | :ref:`mechanics/kane  |
|                       | ``u_dependent =       | <kane_method>` and    |
|                       | [...]``               | the :ref:`kane api.   |
|                       |                       | <kane_lagrange>`      |
|                       | ``u_auxiliary =       |                       |
|                       | [...]``               |                       |
|                       |                       |                       |
|                       | These lists are       |                       |
|                       | passed to the         |                       |
|                       | KanesMethod object.   |                       |
+-----------------------+-----------------------+-----------------------+
| ``Fr()``              | ``KM = KanesMethod(f, | For more details      |
| ``FrStar()``          | q_ind, u_ind, kd_eqs, | refer to              |
|                       | q_dependent, configura| :ref:`mechanics/kane  |
|                       | tion_constraints, u_de| <kane_method>` and    |
|                       | pendent, velocity_cons| the :ref:`kane api.   |
|                       | traints, acceleration_| <kane_lagrange>`      |
|                       | constraints, u_auxilia|                       |
|                       | ry)``                 |                       |
|                       |                       |                       |
|                       | The KanesMethod       |                       |
|                       | object takes a        |                       |
|                       | reference frame       |                       |
|                       | followed by multiple  |                       |
|                       | lists as arguments.   |                       |
|                       |                       |                       |
|                       | ``(fr, frstar) =      |                       |
|                       | KM.kanes_equations(fL,|                       |
|                       | bL)`` where fL and bL |                       |
|                       | are lists of forces   |                       |
|                       | and bodies            |                       |
|                       | respectively.         |                       |
+-----------------------+-----------------------+-----------------------+

Numerical Evaluation and Visualization
----------------------------------------

Autolev’s CODE Option() command allows one to generate Matlab, C, or
Fortran code for numerical evaluation and visualization. Option can be
Dynamics, ODE, Nonlinear or Algebraic.

Numerical evaluation for dynamics can be achieved using PyDy. One can
pass in the KanesMethod object to the System class along with the values
for the constants, specifieds, initial conditions and time steps. The
equations of motion can then be integrated. The plotting is achieved
using matlplotlib. Here is an example from the `PyDy Documentation <https://www.pydy.org/documentation.html>`_
on how it is done::

    from numpy import array, linspace, sin
    from pydy.system import System

    sys = System(kane,
                 constants = {mass: 1.0, stiffness: 1.0,
                              damping: 0.2, gravity: 9.8},
                 specifieds = {force: lambda x, t: sin(t)},
                 initial_conditions = {position: 0.1, speed:-1.0},
                 times = linspace(0.0, 10.0, 1000))

    y = sys.integrate()

    import matplotlib.pyplot as plt
    plt.plot(sys.times, y)
    plt.legend((str(position), str(speed)))
    plt.show()

For information on all the things PyDy can accomplish refer to the
`PyDy Documentation <https://www.pydy.org/documentation.html>`_.

The tools in the PyDy workflow are :

-  `SymPy <https://sympy.org>`_: SymPy is a Python library for
    symbolic computation. It provides computer algebra capabilities
    either as a standalone application, as a library to other
    applications, or live on the web as SymPy Live or SymPy Gamma.

-  `NumPy <https://numpy.org/>`_: NumPy is a library for the
    Python programming language, adding support for large,
    multi-dimensional arrays and matrices, along with a large
    collection of high-level mathematical functions to operate on
    these arrays.

-  `SciPy <https://scipy.org/>`_: SciPy is an open source
    Python library used for scientific computing and technical
    computing. SciPy contains modules for optimization, linear
    algebra, integration, interpolation, special functions, FFT,
    signal and image processing, ODE solvers and other tasks common
    in science and engineering.

-  `IPython <https://ipython.org/>`_: IPython is a command shell
    for interactive computing in multiple programming languages,
    originally developed for the Python programming language, that
    offers introspection, rich media, shell syntax, tab completion,
    and history.

-  `Aesara <https://aesara.readthedocs.io/en/latest/>`_: Aesara is
    a numerical computation library for Python. In Aesara,
    computations are expressed using a NumPy-esque syntax and
    compiled to run efficiently on either CPU or GPU architectures.

-  `Cython <https://cython.org/>`_: Cython is a superset of the
    Python programming language, designed to give C-like performance
    with code that is mostly written in Python. Cython is a compiled
    language that generates CPython extension modules.

-  `matplotlib <https://matplotlib.org/>`_: matplotlib is a
    plotting library for the Python programming language and its
    numerical mathematics extension NumPy.

One will be able to write code equivalent to the Matlab, C or Fortran
code generated by Autolev using these scientific computing tools. It is
recommended to go over these modules to gain an understanding of
scientific computing with Python.

Links
----------

:ref:`SymPy Introductory Tutorial <intro-tutorial>`

:ref:`SymPy Documentation <documentation>`

:ref:`SymPy Physics Vector
Documentation <physics_vector>`

:ref:`SymPy Mechanics
Documentation <physics_mechanics>`

`PyDy Documentation <https://www.pydy.org/documentation.html>`_

`MultiBody Dynamics with Python <https://www.youtube.com/watch?v=mdo2NYtA-xY>`_
