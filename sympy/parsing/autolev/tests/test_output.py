from .... import sin, cos
from .. import parse_autolev


def test_output_01():
    """Autolev example calculates the position, velocity, and accleration of a
    point and expresses in a single reference frame::

          (1) FRAMES C,D,F
          (2) VARIABLES FD'',DC''
          (3) CONSTANTS R,L
          (4) POINTS O,E
          (5) SIMPROT(F,D,1,FD)
       -> (6) F_D = [1, 0, 0; 0, COS(FD), -SIN(FD); 0, SIN(FD), COS(FD)]
          (7) SIMPROT(D,C,2,DC)
       -> (8) D_C = [COS(DC), 0, SIN(DC); 0, 1, 0; -SIN(DC), 0, COS(DC)]
          (9) W_C_F> = EXPRESS(W_C_F>, F)
       -> (10) W_C_F> = FD'*F1> + COS(FD)*DC'*F2> + SIN(FD)*DC'*F3>
          (11) P_O_E>=R*D2>-L*C1>
          (12) P_O_E>=EXPRESS(P_O_E>, D)
       -> (13) P_O_E> = -L*COS(DC)*D1> + R*D2> + L*SIN(DC)*D3>
          (14) V_E_F>=EXPRESS(DT(P_O_E>,F),D)
       -> (15) V_E_F> = L*SIN(DC)*DC'*D1> - L*SIN(DC)*FD'*D2> + (R*FD'+L*COS(DC)*DC')*D3>
          (16) A_E_F>=EXPRESS(DT(V_E_F>,F),D)
       -> (17) A_E_F> = L*(COS(DC)*DC'^2+SIN(DC)*DC'')*D1> + (-R*FD'^2-2*L*COS(DC)*DC'*FD'-L*SIN(DC)*FD'')*D2> + (R*FD''+L*COS(DC)*DC''-L*SIN(DC)*DC'^2-L*SIN(DC)*FD'^2)*D3>

    """
    autolev_input = """\
FRAMES C,D,F
VARIABLES FD'',DC''
CONSTANTS R,L
POINTS O,E
SIMPROT(F,D,1,FD)
SIMPROT(D,C,2,DC)
W_C_F>=EXPRESS(W_C_F>,F)
P_O_E>=R*D2>-L*C1>
P_O_E>=EXPRESS(P_O_E>,D)
V_E_F>=EXPRESS(DT(P_O_E>,F),D)
A_E_F>=EXPRESS(DT(V_E_F>,F),D)\
"""

    sympy_input = parse_autolev(autolev_input)

    g = {}
    l = {}
    exec(sympy_input, g, l)

    w_c_f = l['frame_c'].ang_vel_in(l['frame_f'])
    # P_O_E> means "the position of point E wrt to point O"
    p_o_e = l['point_e'].pos_from(l['point_o'])
    v_e_f = l['point_e'].vel(l['frame_f'])
    a_e_f = l['point_e'].acc(l['frame_f'])

    # NOTE : The Autolev outputs above were manually transformed into
    # equivalent SymPy physics vector expressions. Would be nice to automate
    # this transformation.
    expected_w_c_f = (l['fd'].diff()*l['frame_f'].x +
                      cos(l['fd'])*l['dc'].diff()*l['frame_f'].y +
                      sin(l['fd'])*l['dc'].diff()*l['frame_f'].z)

    assert (w_c_f - expected_w_c_f).simplify() == 0

    expected_p_o_e = (-l['l']*cos(l['dc'])*l['frame_d'].x +
                      l['r']*l['frame_d'].y +
                      l['l']*sin(l['dc'])*l['frame_d'].z)

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
