from .... import sin, cos
from .. import parse_autolev


def test_output_01():
    """
    Autolev example calculates the position, velocity, and accleration of a
    point and expresses in a single reference frame::

          (2) FRAMES C,D,F
          (3) VARIABLES FD'',DC''
          (4) CONSTANTS R,L
          (5) POINTS O,E
          (6) SIMPROT(F,D,1,FD)
       -> (7) F_D = [1, 0, 0; 0, COS(FD), -SIN(FD); 0, SIN(FD), COS(FD)]
          (8) SIMPROT(D,C,2,DC)
       -> (9) D_C = [COS(DC), 0, SIN(DC); 0, 1, 0; -SIN(DC), 0, COS(DC)]
          (10) P_O_E>=R*D2>-L*C1>
       -> (11) P_O_E> = -L*C1> + R*D2>
          (12) V_E_F>=EXPRESS(DT(P_O_E>,F),D)
       -> (13) V_E_F> = L*SIN(DC)*DC'*D1> - L*SIN(DC)*FD'*D2> + (R*FD'+L*COS(DC)*DC')*D3>
          (14) A_E_F>=EXPRESS(DT(V_E_F>,F),D)
       -> (15) A_E_F> = L*(COS(DC)*DC'^2+SIN(DC)*DC'')*D1> + (-R*FD'^2-2*L*COS(DC)*DC'*FD'-L*SIN(DC)*FD'')*D2> + (R*FD''+L*COS(DC)*DC''-L*SIN(DC)*DC'^2-L*SIN(DC)*FD'^2)*D3>
          (16) EXPRESS(P_O_E>, D)
       -> (17) P_O_E> = -L*COS(DC)*D1> + R*D2> + L*SIN(DC)*D3>

    """
    autolev_input = """\
FRAMES C,D,F
VARIABLES FD'',DC''
CONSTANTS R,L
POINTS O,E
SIMPROT(F,D,1,FD)
SIMPROT(D,C,2,DC)
P_O_E>=R*D2>-L*C1>
P_O_E>=EXPRESS(P_O_E>,D)
V_E_F>=EXPRESS(DT(P_O_E>,F),D)
A_E_F>=EXPRESS(DT(V_E_F>,F),D)\
"""

    sympy_input = parse_autolev(autolev_input)

    g = {}
    l = {}
    exec(sympy_input, g, l)

    p_o_e = l['point_o'].pos_from(l['point_e'])
    v_e_f = l['point_e'].vel(l['frame_f'])
    a_e_f = l['point_e'].acc(l['frame_f'])

    # NOTE : The Autolev outputs above were manually transformed into
    # equivalent SymPy physics vector expressions. Would be nice to automate
    # this transformation.
    expected_p_o_e = (-l['l']*cos(l['dc'])*l['frame_d'].z +
                      l['r']*l['frame_d'].y +
                      l['l']*sin(l['dc'])*l['frame_d'].z)

    # TODO : This test fails because the parser seems to output incorrect code.
    # P_O_E> means "the position of point O wrt to point E" (I think).
    assert (p_o_e - expected_p_o_e).simplify() == 0

    expected_v_e_f = (l['l']*sin(l['dc'])*l['dc'].diff()*l['frame_d'].x -
                      l['l']*sin(l['dc'])*l['fd'].diff()*l['frame_d'].y +
                      (l['r']*l['fd'].diff() +
                       l['l']*cos(l['dc'])*l['dc'].diff())*l['frame_d'].z)
    assert (v_e_f - expected_v_e_f).simplify() == 0

    expected_a_e_f = (l['l']*(cos(l['dc'])*l['dc'].diff()**2 +
                              sin(l['dc'])*l['dc'].diff().diff())*l['frame_d'].x +
                      (-l['r']*l['fd'].diff()**2 -
                       2*l['l']*cos(l['dc'])*l['dc'].diff()*l['fd'].diff() -
                       l['l']*sin(l['dc'])*l['fd'].diff().diff())*l['frame_d'].y +
                      (l['r']*l['fd'].diff().diff() +
                       l['l']*cos(l['dc'])*l['dc'].diff().diff() -
                       l['l']*sin(l['dc'])*l['dc'].diff()**2 -
                       l['l']*sin(l['dc'])*l['fd'].diff()**2)*l['frame_d'].z)
    assert (a_e_f - expected_a_e_f).simplify() == 0
