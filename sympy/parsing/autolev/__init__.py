from sympy.external import import_module


def parse_autolev(input, output=None, include_pydy=False):
    """Parses Autolev code (version 4.1) to SymPy code.

    Parameters
    ----------
    input: str
         1. Can be the name of a file containing Autolev code.
         2. Can be a string containing Autolev code.

    output: str
          1. Can be the name of an output file to which the SymPy code should be written to.
          2. Can be the string "list". In this it returns a list containing the SymPy code.
             Each element in the list corresponds to one line of code.
          3. If nothing is passed in the parsed code is printed out to stdout.

    include_pydy: boolean
                The parser will output PyDy ode integration code (if the Autolev code calls for it)
                if this is set to True. It is set to False by default.


    Example (Double Pendulum)
    -------
    Suppose this is the file we need to parse:
    ------------FILE--------------
    % double_pendulum.al
    MOTIONVARIABLES' Q{2}', U{2}'
    CONSTANTS L,M,G
    NEWTONIAN N
    FRAMES A,B
    SIMPROT(N, A, 3, Q1)
    SIMPROT(N, B, 3, Q2)
    W_A_N>=U1*N3>
    W_B_N>=U2*N3>
    POINT O
    PARTICLES P,R
    P_O_P> = L*A1>
    P_P_R> = L*B1>
    V_O_N> = 0>
    V2PTS(N, A, O, P)
    V2PTS(N, B, P, R)
    MASS P=M, R=M
    Q1' = U1
    Q2' = U2
    GRAVITY(G*N1>)
    ZERO = FR() + FRSTAR()
    KANE()
    INPUT M=1,G=9.81,L=1
    INPUT Q1=.1,Q2=.2,U1=0,U2=0
    INPUT TFINAL=10, INTEGSTP=.01
    CODE DYNAMICS() some_filename.c
    ------------FILE--------------
    Here parse_autolev() is used without passing in a parameter for the output file
    so it defaults to stdout.
    >>> from sympy.parsing.autolev import parse_autolev  # doctest: +SKIP
    >>> parse_autolev("double_pendulum.al", include_pydy=True)  # doctest: +SKIP
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
    force_p = particle_p.mass*(g*frame_n.x)
    force_r = particle_r.mass*(g*frame_n.x)
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

    y=sys.integrate()
    """

    _autolev = import_module(
        'sympy.parsing.autolev._parse_autolev_antlr',
        __import__kwargs={'fromlist': ['X']})

    if _autolev is not None:
        return _autolev.parse_autolev(input, output, include_pydy)
