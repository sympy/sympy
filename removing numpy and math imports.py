import sympy.physics.mechanics as me
import sympy as sm

frame_c = me.ReferenceFrame('c')
frame_d = me.ReferenceFrame('d')
frame_f = me.ReferenceFrame('f')
fd, dc = me.dynamicsymbols('fd dc')
fdd, dcd = me.dynamicsymbols('fd dc', 1)
fdd2, dcd2 = me.dynamicsymbols('fd dc', 2)
r, l = sm.symbols('r l', real=True)
point_o = me.Point('o')
point_e = me.Point('e')
frame_d.orient(frame_f, 'Axis', [fd, frame_f.x])
frame_c.orient(frame_d, 'Axis', [dc, frame_d.y])
frame_c.set_ang_vel(frame_f, (frame_c.ang_vel_in(frame_f)).express(frame_f))
point_e.set_pos(point_o, r*frame_d.y-l*frame_c.x)
point_e.set_pos(point_o, (point_e.pos_from(point_o)).express(frame_d))
point_e.set_vel(frame_f, ((point_e.pos_from(point_o)).dt(frame_f)).express(frame_d))
point_e.set_acc(frame_f, ((point_e.vel(frame_f)).dt(frame_f)).express(frame_d))
