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

The parser has been built using the `ANTLR <http://www.antlr.org/>`_ framework and its main purpose
is to help former users of Autolev to get familiarized with multibody dynamics
in SymPy.

The sections below shall discuss details of the parser like usage, gotchas,
issues and future improvements.
For a detailed comparison of Autolev and SymPy/PyDy you might want to look at
the PyDy for Autolev Users guide.

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
    >>> parse_autolev('double_pendulum.al', 'double_pendulum.py', include_pydy=True)
    # The second parameter (to specify the output file) is optional. The
    # default action is to print the output code to stdout.
    #
    # The include_pydy flag is False by default. Setting it to True will
    # allow the PyDy code to be outputed if applicable.

This is the parsed code in ``double_pendulum.py``::

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
If you are completely new to SymPy mechanics, the PyDy for Autolev Users
guide should help. You might also have to use basic SymPy simplifications
and manipulations like ``trigsimp()``, ``expand()``, ``evalf()`` etc for 
getting outputs similar to Autolev.
Refer to the `SymPy Tutorial <http://docs.sympy.org/latest/tutorial/index.html>`_
to know more about these.
