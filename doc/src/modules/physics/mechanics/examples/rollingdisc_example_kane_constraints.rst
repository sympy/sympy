========================================================
A rolling disc, with Kane's method and constraint forces
========================================================

We will now revisit the rolling disc example, except this time we are bringing
the non-contributing (constraint) forces into evidence. See [Kane1985]_ for a
more thorough explanation of this. Here, we will turn on the automatic
simplifcation done when doing vector operations. It makes the outputs nicer for
small problems, but can cause larger vector operations to hang. ::

  >>> from sympy import symbols, sin, cos, tan
  >>> from sympy.physics.mechanics import *
  >>> mechanics_printing(pretty_print=False)
  >>> q1, q2, q3, u1, u2, u3  = dynamicsymbols('q1 q2 q3 u1 u2 u3')
  >>> q1d, q2d, q3d, u1d, u2d, u3d = dynamicsymbols('q1 q2 q3 u1 u2 u3', 1)
  >>> r, m, g = symbols('r m g')

These two lines introduce the extra quantities needed to find the constraint
forces. ::

  >>> u4, u5, u6, f1, f2, f3 = dynamicsymbols('u4 u5 u6 f1 f2 f3')

Most of the main code is the same as before. ::

  >>> N = ReferenceFrame('N')
  >>> Y = N.orientnew('Y', 'Axis', [q1, N.z])
  >>> L = Y.orientnew('L', 'Axis', [q2, Y.x])
  >>> R = L.orientnew('R', 'Axis', [q3, L.y])
  >>> w_R_N_qd = R.ang_vel_in(N)
  >>> R.set_ang_vel(N, u1 * L.x + u2 * L.y + u3 * L.z)

The definition of rolling without slip necessitates that the velocity of the
contact point is zero; as part of bringing the constraint forces into evidence,
we have to introduce speeds at this point, which will by definition always be
zero. They are normal to the ground, along the path which the disc is rolling,
and along the ground in an perpendicular direction. ::

  >>> C = Point('C')
  >>> C.set_vel(N, u4 * L.x + u5 * (Y.z ^ L.x) + u6 * Y.z)
  >>> Dmc = C.locatenew('Dmc', r * L.z)
  >>> vel = Dmc.v2pt_theory(C, N, R)
  >>> I = inertia(L, m / 4 * r**2, m / 2 * r**2, m / 4 * r**2)
  >>> kd = [dot(R.ang_vel_in(N) - w_R_N_qd, uv) for uv in L]

Just as we previously introduced three speeds as part of this process, we also
introduce three forces; they are in the same direction as the speeds, and
represent the constraint forces in those directions. ::

  >>> ForceList = [(Dmc, - m * g * Y.z), (C, f1 * L.x + f2 * (Y.z ^ L.x) + f3 * Y.z)]
  >>> BodyD = RigidBody('BodyD', Dmc, R, m, (I, Dmc))
  >>> BodyList = [BodyD]

  >>> KM = KanesMethod(N, q_ind=[q1, q2, q3], u_ind=[u1, u2, u3], kd_eqs=kd,
  ...           u_auxiliary=[u4, u5, u6])
  >>> (fr, frstar) = KM.kanes_equations(BodyList, ForceList)
  >>> MM = KM.mass_matrix
  >>> forcing = KM.forcing
  >>> rhs = MM.inv() * forcing
  >>> kdd = KM.kindiffdict()
  >>> rhs = rhs.subs(kdd)
  >>> rhs.simplify()
  >>> mprint(rhs)
  Matrix([
  [(4*g*sin(q2) + 6*r*u2*u3 - r*u3**2*tan(q2))/(5*r)],
  [                                       -2*u1*u3/3],
  [                          (-2*u2 + u3*tan(q2))*u1]])
  >>> from sympy import trigsimp, signsimp, collect, factor_terms
  >>> def simplify_auxiliary_eqs(w):
  ...     return signsimp(trigsimp(collect(collect(factor_terms(w), f2), m*r)))
  >>> mprint(KM.auxiliary_eqs.applyfunc(simplify_auxiliary_eqs))
  Matrix([
  [                                      -m*r*(u1*u3 + u2') + f1],
  [-m*r*u1**2*sin(q2) - m*r*u2*u3/cos(q2) + m*r*cos(q2)*u1' + f2],
  [                -g*m + m*r*(u1**2*cos(q2) + sin(q2)*u1') + f3]])
