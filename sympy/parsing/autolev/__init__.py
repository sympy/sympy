from sympy.external import import_module


def parse_autolev(inp, output=None):
    """Parses Autolev code to SymPy code.

    Parameters
    ----------
    inp: str
         1. Can be the name of a file containing Autolev code.
         2. Can be a string containing Autolev code.

    output: str
            1. Can be the name of an output file to which the SymPy code should be written to.
            2. Can be the string "print". In this case the SymPy code is written to stdout.
            3. Can be the string "list". In this it returns a list containing the SymPy code.
                Each element in the list corresponds to one line of code.

    Example (Double Pendulum)
    -------
    >>> from sympy.parsing.autolev import parse_autolev  # doctest: +SKIP
    >>> l = []  # doctest: +SKIP
    >>> l.append("MOTIONVARIABLES' Q{2}', U{2}')"  # doctest: +SKIP
    >>> l.append("CONSTANTS L,M,G)"  # doctest: +SKIP
    >>> l.append("NEWTONIAN N)"  # doctest: +SKIP
    >>> l.append("FRAMES A,B)"  # doctest: +SKIP
    >>> l.append("SIMPROT(N, A, 3, Q1))"  # doctest: +SKIP
    >>> l.append("SIMPROT(N, B, 3, Q2))"  # doctest: +SKIP
    >>> l.append("W_A_N>=U1*N3>)"  # doctest: +SKIP
    >>> l.append("W_B_N>=U2*N3>)"  # doctest: +SKIP
    >>> l.append("POINT O)"  # doctest: +SKIP
    >>> l.append("PARTICLES P,R)"  # doctest: +SKIP
    >>> l.append("P_O_P> = L*A1>)"  # doctest: +SKIP
    >>> l.append("P_P_R> = L*B1>)"  # doctest: +SKIP
    >>> l.append("V_O_N> = 0>)"  # doctest: +SKIP
    >>> l.append("V2PTS(N, A, O, P))"  # doctest: +SKIP
    >>> l.append("V2PTS(N, B, P, R))"  # doctest: +SKIP
    >>> l.append("MASS P=M, R=M)"  # doctest: +SKIP
    >>> l.append("Q1' = U1)"  # doctest: +SKIP
    >>> l.append("Q2' = U2)"  # doctest: +SKIP
    >>> l.append("GRAVITY(G*N1>))"  # doctest: +SKIP
    >>> l.append("ZERO = FR() + FRSTAR())"  # doctest: +SKIP
    >>> l.append("KANE())"  # doctest: +SKIP
    >>> l.append("INPUT M=1,G=9.81,L=1)"  # doctest: +SKIP
    >>> l.append("INPUT Q1=.1,Q2=.2,U1=0,U2=0)"  # doctest: +SKIP
    >>> l.append("INPUT TFINAL=10, INTEGSTP=.01)"  # doctest: +SKIP
    >>> l.append("CODE DYNAMICS() double_pendulum.c)"  # doctest: +SKIP
    >>> parse_autolev("\\n".join(l), "print")  # doctest: +SKIP
    import sympy.physics.mechanics as me
    import sympy as sm
    import math as m
    import numpy as np
    <BLANKLINE>
    q1, q2, u1, u2 = me.dynamicsymbols('q1 q2 u1 u2')
    q1d, q2d, u1d, u2d = me.dynamicsymbols('q1 q2 u1 u2', 1)
    l, m, g=sm.symbols('l m g', real=True)
    frame_n=me.ReferenceFrame('n')
    frame_a=me.ReferenceFrame('a')
    frame_b=me.ReferenceFrame('b')
    frame_a.orient(frame_n, 'Axis', [q1, frame_n.z])
    frame_b.orient(frame_n, 'Axis', [q2, frame_n.z])
    frame_a.set_ang_vel(frame_n, u1*frame_n.z)
    frame_b.set_ang_vel(frame_n, u2*frame_n.z)
    point_o=me.Point('o')
    particle_p=me.Particle('p', me.Point('p_pt'), sm.Symbol('m'))
    particle_r=me.Particle('r', me.Point('r_pt'), sm.Symbol('m'))
    particle_p.point.set_pos(point_o, l*frame_a.x)
    particle_r.point.set_pos(particle_p.point, l*frame_b.x)
    point_o.set_vel(frame_n, 0)
    particle_p.point.v2pt_theory(point_o,frame_n,frame_a)
    particle_r.point.v2pt_theory(particle_p.point,frame_n,frame_b)
    particle_p.mass = m
    particle_r.mass = m
    kd_eqs = [q1d - u1, q2d - u2]
    forceList = [(particle_p.point,particle_p.mass*(g*frame_n.x)), (particle_r.point,particle_r.mass*(g*frame_n.x))]
    kane = me.KanesMethod(frame_n, q_ind=[q1,q2], u_ind=[u1, u2], kd_eqs = kd_eqs)
    fr, frstar = kane.kanes_equations([particle_p, particle_r], forceList)
    zero = fr+frstar
    from pydy.system import System
    sys = System(kane, constants = {l:1, m:1, g:9.81},
    specifieds={},
    initial_conditions={q1:.1, q2:.2, u1:0, u2:0},
    times = np.linspace(0.0, 10, 10/.01))
    <BLANKLINE>
    y=sys.integrate()
    <BLANKLINE>
    """

    _autolev = import_module(
        'sympy.parsing.autolev._parse_autolev_antlr',
        __import__kwargs={'fromlist': ['X']})

    if _autolev is not None:
        return _autolev.parse_autolev(inp, output)
