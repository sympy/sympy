.. _autolev_parser:

==============
Autolev Parser
==============

.. role:: input(strong)

Introduction
============
Autolev (now superseded by MotionGenesis) is a domain specific language
used for symbolic multibody dynamics. The SymPy mechanics module now has
enough power and functionality to be a fully featured symbolic dynamics
module. This parser parses Autolev (version 4.1) code to SymPy code by making
use of SymPy's math libraries and the mechanics module.

The parser has been built using the `ANTLR <https://www.antlr.org/>`_ framework and its main purpose
is to help former users of Autolev to get familiarized with multibody dynamics
in SymPy.

The sections below shall discuss details of the parser like usage, gotchas,
issues and future improvements.
For a detailed comparison of Autolev and SymPy Mechanics you might want to look at
the :ref:`SymPy Mechanics for Autolev Users guide <sympy_mechanics_for_autolev_users>`.

.. _usage:

Usage
=====

We first start with an Autolev code file.

Let us take this example
(Comments ``%`` have been included to show the Autolev responses):

.. code-block:: none

   % double_pendulum.al
   %-------------------
   MOTIONVARIABLES' Q{2}', U{2}'
   CONSTANTS L,M,G
   NEWTONIAN N
   FRAMES A,B
   SIMPROT(N, A, 3, Q1)
   % -> N_A = [COS(Q1), -SIN(Q1), 0; SIN(Q1), COS(Q1), 0; 0, 0, 1]
   SIMPROT(N, B, 3, Q2)
   % -> N_B = [COS(Q2), -SIN(Q2), 0; SIN(Q2), COS(Q2), 0; 0, 0, 1]
   W_A_N>=U1*N3>
   % -> W_A_N> = U1*N3>
   W_B_N>=U2*N3>
   % -> W_B_N> = U2*N3>
   POINT O
   PARTICLES P,R
   P_O_P> = L*A1>
   % -> P_O_P> = L*A1>
   P_P_R> = L*B1>
   % -> P_P_R> = L*B1>
   V_O_N> = 0>
   % -> V_O_N> = 0>
   V2PTS(N, A, O, P)
   % -> V_P_N> = L*U1*A2>
   V2PTS(N, B, P, R)
   % -> V_R_N> = L*U1*A2> + L*U2*B2>
   MASS P=M, R=M
   Q1' = U1
   Q2' = U2
   GRAVITY(G*N1>)
   % -> FORCE_P> = G*M*N1>
   % -> FORCE_R> = G*M*N1>
   ZERO = FR() + FRSTAR()
   % -> ZERO[1] = -L*M*(2*G*SIN(Q1)+L*(U2^2*SIN(Q1-Q2)+2*U1'+COS(Q1-Q2)*U2'))
   % -> ZERO[2] = -L*M*(G*SIN(Q2)-L*(U1^2*SIN(Q1-Q2)-U2'-COS(Q1-Q2)*U1'))
   KANE()
   INPUT M=1,G=9.81,L=1
   INPUT Q1=.1,Q2=.2,U1=0,U2=0
   INPUT TFINAL=10, INTEGSTP=.01
   CODE DYNAMICS() some_filename.c

The parser can be used as follows::

    >>> from sympy.parsing.autolev import parse_autolev
    >>> sympy_code = parse_autolev(open('double_pendulum.al'), include_numeric=True)

    # The include_pydy flag is False by default. Setting it to True will
    # enable PyDy simulation code to be outputted if applicable.

    >>> print(sympy_code)
    import sympy.physics.mechanics as me
    import sympy as sm
    import math as m
    import numpy as np

    q1, q2, u1, u2 = me.dynamicsymbols('q1 q2 u1 u2')
    q1d, q2d, u1d, u2d = me.dynamicsymbols('q1 q2 u1 u2', 1)
    l, m, g=sm.symbols('l m g', real=True)
    frame_n=me.ReferenceFrame('n')
    frame_a=me.ReferenceFrame('a')
    frame_b=me.ReferenceFrame('b')
    frame_a.orient(frame_n, 'Axis', [q1, frame_n.z])
    # print(frame_n.dcm(frame_a))
    frame_b.orient(frame_n, 'Axis', [q2, frame_n.z])
    # print(frame_n.dcm(frame_b))
    frame_a.set_ang_vel(frame_n, u1*frame_n.z)
    # print(frame_a.ang_vel_in(frame_n))
    frame_b.set_ang_vel(frame_n, u2*frame_n.z)
    # print(frame_b.ang_vel_in(frame_n))
    point_o=me.Point('o')
    particle_p=me.Particle('p', me.Point('p_pt'), sm.Symbol('m'))
    particle_r=me.Particle('r', me.Point('r_pt'), sm.Symbol('m'))
    particle_p.point.set_pos(point_o, l*frame_a.x)
    # print(particle_p.point.pos_from(point_o))
    particle_r.point.set_pos(particle_p.point, l*frame_b.x)
    # print(particle_p.point.pos_from(particle_r.point))
    point_o.set_vel(frame_n, 0)
    # print(point_o.vel(frame_n))
    particle_p.point.v2pt_theory(point_o,frame_n,frame_a)
    # print(particle_p.point.vel(frame_n))
    particle_r.point.v2pt_theory(particle_p.point,frame_n,frame_b)
    # print(particle_r.point.vel(frame_n))
    particle_p.mass = m
    particle_r.mass = m
    force_p = particle_p.mass*(g*frame_n.x)
    # print(force_p)
    force_r = particle_r.mass*(g*frame_n.x)
    # print(force_r)
    kd_eqs = [q1d - u1, q2d - u2]
    forceList = [(particle_p.point,particle_p.mass*(g*frame_n.x)), (particle_r.point,particle_r.mass*(g*frame_n.x))]
    kane = me.KanesMethod(frame_n, q_ind=[q1,q2], u_ind=[u1, u2], kd_eqs = kd_eqs)
    fr, frstar = kane.kanes_equations([particle_p, particle_r], forceList)
    zero = fr+frstar
    # print(zero)
    #---------PyDy code for integration----------
    from pydy.system import System
    sys = System(kane, constants = {l:1, m:1, g:9.81},
    specifieds={},
    initial_conditions={q1:.1, q2:.2, u1:0, u2:0},
    times = np.linspace(0.0, 10, 10/.01))

    y=sys.integrate()


The commented code is not part of the output code. The print
statements demonstrate how to get responses similar to the ones in the
Autolev file.
Note that we need to use SymPy functions like ``.ang_vel_in()``, ``.dcm()``
etc in many cases unlike directly printing out the variables like ``zero``.
If you are completely new to SymPy mechanics, the :ref:`SymPy Mechanics for Autolev Users guide <sympy_mechanics_for_autolev_users>`
guide should help. You might also have to use basic SymPy simplifications
and manipulations like ``trigsimp()``, ``expand()``, ``evalf()`` etc for
getting outputs similar to Autolev.
Refer to the `SymPy Tutorial <https://docs.sympy.org/latest/tutorial/index.html>`_
to know more about these.

.. _gotchas_autolev:

Gotchas
=======

- Don't use variable names that conflict with Python's reserved words.
  This is one example where this is violated:

  .. code-block:: none

     %Autolev Code
     %------------
     LAMBDA = EIG(M)

  .. code-block:: python

     #SymPy Code
     #----------
     lambda = sm.Matrix([i.evalf() for i in (m).eigenvals().keys()])

------------------------------------------------------------------------

- Make sure that the names of vectors and scalars are different.
  Autolev treats these differently but these will get overwritten in Python.
  The parser currently allows the names of bodies and scalars/vectors to
  coincide but doesn't do this between scalars and vectors.
  This should probably be changed in the future.

  .. code-block:: none

     %Autolev Code
     %------------
     VARIABLES X,Y
     FRAMES A
     A> = X*A1> + Y*A2>
     A = X+Y

  .. code-block:: python

     #SymPy Code
     #----------
     x, y = me.dynamicsymbols('x y')
     frame_a = me.ReferenceFrame('a')
     a = x*frame_a.x + y*frame_a.y
     a = x + y
     # Note how frame_a is named differently so it doesn't cause a problem.
     # On the other hand, 'a' gets rewritten from a scalar to a vector.
     # This should be changed in the future.


------------------------------------------------------------------------

- When dealing with Matrices returned by functions, one must check the
  order of the values as they may not be the same as in Autolev. This is
  especially the case for eigenvalues and eigenvectors.

  .. code-block:: none

     %Autolev Code
     %------------
     EIG(M, E1, E2)
     % -> [5; 14; 13]
     E2ROW = ROWS(E2, 1)
     EIGVEC> = VECTOR(A, E2ROW)

  .. code-block:: python

     #SymPy Code
     #----------
     e1 = sm.Matrix([i.evalf() for i in m.eigenvals().keys()])
     # sm.Matrix([5;13;14]) different order
     e2 = sm.Matrix([i[2][0].evalf() for i in m.eigenvects()]).reshape(m.shape[0], m.shape[1])
     e2row = e2.row(0)
     # This result depends on the order of the vectors in the eigenvecs.
     eigenvec = e2row[0]*a.x + e2row[1]*a.y + e2row[2]*a.y

------------------------------------------------------------------------

- When using ``EVALUATE``, use something like ``90*UNITS(deg,rad)`` for
  angle substitutions as radians are the default in SymPy.
  You could also add ``np.deg2rad()`` directly in the SymPy code.

  This need not be done for the output code (generated on parsing the
  ``CODE`` commands) as the parser takes care of this when ``deg`` units
  are given in the ``INPUT`` declarations.

  The ``DEGREES`` setting, on the other hand, works only in some cases like
  in ``SIMPROT`` where an angle is expected.

  .. code-block:: none

     %Autolev Code
     %------------
     A> = Q1*A1> + Q2*A2>
     B> = EVALUATE(A>, Q1:30*UNITS(DEG,RAD))

  .. code-block:: python

     #SymPy Code
     #----------
     a = q1*a.frame_a.x + q2*frame_a.y
     b = a.subs({q1:30*0.0174533})
     # b = a.subs({q1:np.deg2rad(30)}

------------------------------------------------------------------------

- Most of the Autolev settings have not been parsed and have no effect on the parser.
  The only ones that work somewhat are ``COMPLEX`` and ``DEGREES``.
  It is advised to look into alternatives to these in SymPy and Python.

------------------------------------------------------------------------

- The ``REPRESENT`` command is not supported.
  Use the ``MATRIX``, ``VECTOR`` or ``DYADIC`` commands instead.
  Autolev 4.1 suggests these over ``REPRESENT`` as well while still allowing
  it but the parser doesn't parse it.

------------------------------------------------------------------------

- Do not use variables declarations of the type ``WO{3}RD{2,4}``.
  The parser can only handle one variable name followed by one pair
  of curly braces and any number of ``'`` s.
  You would have to declare all the cases manually if you want to achieve
  something like ``WO{3}RD{2,4}``.

------------------------------------------------------------------------

- The parser can handle normal versions of most commands but it may not
  parse functions with Matrix arguments properly in most cases.
  Eg:

  ``M=COEF([E1;E2],[U1,U2,U3])``

  This would compute the coefficients of ``U1``, ``U2`` and ``U3`` in ``E1``
  and ``E2``. It is preferable to manually construct a Matrix using the
  regular versions of these commands.

  .. code-block:: none

     %Autolev Code
     %------------
     % COEF([E1;E2],[U1,U2,U3])
     M = [COEF(E1,U1),COEF(E1,U2),COEF(E1,U3) &
         ;COEF(E2,U1),COEF(E2,U2),COEF(E2,U3)]

------------------------------------------------------------------------

- ``MOTIONVARIABLE`` declarations must be used for the generalized coordinates
  and speeds and all other variables must be declared in regular
  ``VARIABLE`` declarations.
  The parser requires this to distinguish between them to pass the correct
  parameters to the Kane's method object.

  It is also preferred to always declare the speeds corresponding to the
  coordinates and to pass in the kinematic differential equations.
  The parser is able to handle some cases where this isn't the case by
  introducing some dummy variables of its own but SymPy on its own
  does require them.

  Also note that older Autolev declarations like ``VARIABLES U{3}'`` are not
  supported either.

  .. code-block:: none

     %Autolev Code
     %------------
     MOTIONVARIABLES' Q{2}', U{2}'
     % ----- OTHER LINES ----
     Q1' = U1
     Q2' = U2
     ----- OTHER LINES ----
     ZERO = FR() + FRSTAR()

  .. code-block:: python

     #SymPy Code
     #----------
     q1, q2, u1, u2 = me.dynamicsymbols('q1 q2 u1 u2')
     q1d, q2d, u1d, u2d = me.dynamicsymbols('q1 q2 u1 u2', 1)

     # ------- other lines -------

     kd_eqs = [q1d - u1, q2d - u2]
     kane = me.KanesMethod(frame_n, q_ind=[q1,q2], u_ind=[u1, u2], kd_eqs = kd_eqs)
     fr, frstar = kane.kanes_equations([particle_p, particle_r], forceList)
     zero = fr+frstar

------------------------------------------------------------------------

- Need to change ``me.dynamicsymbols._t`` to ``me.dynamicsymbols('t')`` for
  all occurrences of it in the Kane's equations. For example have a look at
  line 10 of this `spring damper example <https://github.com/sympy/sympy/blob/master/sympy/parsing/autolev/test-examples/pydy-example-repo/mass_spring_damper.py#L10>`_.
  This equation is used in forming the Kane's equations so we need to
  change ``me.dynamicsymbols._t`` to ``me.dynamicsymbols('t')`` in this case.

  The main reason that this needs to be done is because PyDy
  requires time dependent specifieds to be explicitly laid out while
  Autolev simply takes care of the stray time variables in the equations
  by itself.

  The problem is that PyDy's System class does not accept
  ``dynamicsymbols._t`` as a specified. Refer to issue `#396 <https://github.com/pydy/pydy/issues/396>`_.
  This change is not actually ideal so a better solution should be figured
  out in the future.

------------------------------------------------------------------------

- The parser creates SymPy ``symbols`` and ``dynamicsymbols`` by parsing
  variable declarations in the Autolev Code.

  For intermediate expressions which are directly initialized the parser
  does not create SymPy symbols. It just assigns them to the expression.

  On the other hand, when a declared variable is assigned to an expression,
  the parser stores the expression against the variable in a dictionary so
  as to not reassign it to a completely different entity. This constraint
  is due to the inherent nature of Python and how it differs from a language
  like Autolev.

  Also, Autolev seems to be able to assume whether to use a variable or the
  rhs expression that variable has been assigned to in equations even
  without an explicit ``RHS()`` call in some cases.
  For the parser to work correctly however, it is better to use ``RHS()``
  wherever a variable's rhs expression is meant to be used.

  .. code-block:: none

     %Autolev Code
     %------------
     VARIABLES X, Y
     E = X + Y
     X = 2*Y

     RHS_X = RHS(X)

     I1 = X
     I2 = Y
     I3 = X + Y

     INERTIA B,I1,I2,I3
     % -> I_B_BO>> = I1*B1>*B1> + I2*B2>*B2> + I3*B3>*B3>

  .. code-block:: python

     #SymPy Code
     #----------
     x,y = me.dynamicsymbols('x y')
     e = x + y  # No symbol is made out of 'e'

     # an entry like {x:2*y} is stored in an rhs dictionary

     rhs_x = 2*y

     i1 = x  # again these are not made into SymPy symbols
     i2 = y
     i3 = x + y

     body_b.inertia = (me.inertia(body_b_f, i1, i2, i3), b_cm)
     # This prints as:
     # x*b_f.x*b_f.x + y*b_f.y*b_f.y + (x+y)*b_f.z*b_f.z
     # while Autolev's output has I1,I2 and I3 in it.
     # Autolev however seems to know when to use the RHS of I1,I2 and I3
     # based on the context.

------------------------------------------------------------------------

- This is how the ``SOLVE`` command is parsed:

  .. code-block:: none

     %Autolev Code
     %------------
     SOLVE(ZERO,X,Y)
     A = RHS(X)*2 + RHS(Y)

  .. code-block:: python

     #SymPy Code
     #----------
     print(sm.solve(zero,x,y))
     # Behind the scenes the rhs of x
     # is set to sm.solve(zero,x,y)[x].
     a = sm.solve(zero,x,y)[x]*2 + sm.solve(zero,x,y)[y]

  The indexing like ``[x]`` and ``[y]`` doesn't always work so you might want to
  look at the underlying dictionary that solve returns and index it correctly.

------------------------------------------------------------------------

- Inertia declarations and Inertia functions work somewhat differently in
  the context of the parser. This might be hard to understand at first
  but this had to be done to bridge the gap due to the differences in
  SymPy and Autolev. Here are some points about them:

  1. Inertia declarations (``INERTIA B,I1,I2,I3``) set the inertias of rigid
  bodies.

  2. Inertia setters of the form ``I_C_D>> = expr`` however, set the inertias
  only when C is a body. If C is a particle then ``I_C_D>> = expr``
  simply parses to ``i_c_d = expr`` and ``i_c_d`` acts like a regular variable.

  3. When it comes to inertia getters (``I_C_D>>`` used in an expression or
  ``INERTIA`` commands), these MUST be used with the ``EXPRESS`` command
  to specify the frame as SymPy needs this information to compute the
  inertia dyadic.

  .. code-block:: none

     %Autolev Code
     %------------
     INERTIA B,I1,I2,I3
     I_B_BO>> = X*A1>*A1> + Y*A2>*A2>  % Parser will set the inertia of B
     I_P_Q>> = X*A1>*A1> + Y^2*A2>*A2> % Parser just parses it as i_p_q = expr

     E1 = 2*EXPRESS(I_B_O>>,A)
     E2 =  I_P_Q>>
     E3 = EXPRESS(I_P_O>>,A)
     E4 = EXPRESS(INERTIA(O),A)

     % In E1 we are using the EXPRESS command with I_B_O>> which makes
     % the parser and SymPy compute the inertia of Body B about point O.

     % In E2 we are just using the dyadic object I_P_Q>> (as I_P_Q>> = expr
     % doesn't act as a setter) defined above and not asking the parser
     % or SymPy to compute anything.

     % E3 asks the parser to compute the inertia of P about point O.

     % E4 asks the parser to compute the inertias of all bodies wrt about O.

------------------------------------------------------------------------

- In an inertia declaration of a body, if the inertia is being set about
  a point other than the center of mass, one needs to make sure that
  the position vector setter for that point and the center of mass appears
  before the inertia declaration as SymPy will throw an error otherwise.

  .. code-block:: none

     %Autolev Code
     %------------
     P_SO_O> = X*A1>
     INERTIA S_(O) I1,I2,I3

------------------------------------------------------------------------

- Note that all Autolev commands have not been implemented. The parser
  now covers the important ones in their basic forms. If you are doubtful
  whether a command is included or not, please have a look at `this file <https://github.com/sympy/sympy/blob/master/sympy/parsing/autolev/_listener_autolev_antlr.py>`_
  in the source code. Search for "<command>" to verify this. Looking at the
  code for the specific command will also give an idea about what form it
  is expected to work in.

.. _issues:

Limitations and Issues
======================

- A lot of the issues have already been discussed in the Gotchas section.
  Some of these are:

  - Vector names coinciding with scalar names are overwritten in Python.
  - Some convenient variable declarations aren't parsed.
  - Some convenient forms of functions to return matrices aren't parsed.
  - Settings aren't parsed.
  - symbols and rhs expressions work very differently in Python which might
    cause undesirable results.
  - Dictionary indexing for the parsed code of the ``SOLVE`` command is
    not proper in many cases.
  - Need to change ``dynamicsymbols._t`` to ``dynamicsymbols('t')`` for the
    PyDy simulation code to work properly.

Here are some other ones:

- Eigenvectors do not seem to work as expected. The values in Autolev and SymPy
  are not the same in many cases.

- Block matrices aren't parsed by the parser. It would actually be easier
  to make a change in SymPy to allow matrices to accept other matrices for
  arguments.

- The SymPy equivalent of the ``TAYLOR`` command ``.series()`` does not work
  with ``dynamicsymbols()``.

- Only ``DEPENDENT`` constraints are currently parsed. Need to parse
  ``AUXILIARY`` constraints as well. This should be done soon as it isn't
  very difficult.

- None of the energy and momentum functions are parsed right now. It would
  be nice to get these working as well. Some changes should probably be made
  to SymPy. For instance, SymPy doesn't have a function equivalent to ``NICHECK()``.

- The numerical integration parts work properly only in the case of the
  ``KANE`` command with no arguments. Things like ``KANE(F1,F2)`` do not currently
  work.

- Also, the PyDy numerical simulation code works only for cases where the
  matrix say ``ZERO = FR() + FRSTAR()`` is solved for. It doesn't work well when the
  matrix has some other equations plugged in as well. One hurdle
  faced in achieving this was that PyDy's System class automatically takes
  in the ``forcing_full`` and ``mass_matrix_full`` and solves them without giving the
  user the flexibility to specify the equations. It would be nice to add
  this functionality to the System class.


.. _future_improvements:

Future Improvements
===================

1. Completing Dynamics Online
-----------------------------
The parser has been built by referring to and parsing codes from the
`Autolev Tutorial <https://mae.ufl.edu/~fregly/PDFs/autolev_tutorial.pdf>`_
and the book *Dynamics Online: Theory and Implementation Using Autolev*.
Basically, the process involved going through each of these codes,
validating the parser results and improving the rules if required
to make sure the codes parsed well.

The parsed codes of these are available on GitLab `here <https://gitlab.com/sympy/autolev-test-examples>`_.
The repo is private so access needs to be requested.
As of now, most codes till Chapter 4 of *Dynamics Online* have been parsed.

Completing all the remaining codes of the book (namely, *2-10*, *2-11*, *rest
of Ch4*, *Ch5* and *Ch6* (less important) ) would make the parser more complete.


2. Fixing Issues
----------------
The second thing to do would be to go about fixing the problems described
above in the :ref:`Gotchas <gotchas_autolev>` and :ref:`Limitations and Issues <issues>`
sections in order of priority and ease. Many of these require changes
in the parser code while some of these are better fixed by adding some
functionality to SymPy.


3. Switching to an AST
----------------------
The parser is currently built using a kind of Concrete Syntax Tree (CST)
using the `ANTLR <https://www.antlr.org/>`_ framework. It would be ideal to switch from a CST to an
Abstract Syntax Tree (AST). This way, the parser code will be independent
of the ANTLR grammar which makes it a lot more flexible. It would also be
easier to make changes to the grammar and the rules of the parser.
