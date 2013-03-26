
from numpy import mgrid,shape,swapaxes,zeros,log,exp,sin,cos,tan
 
eps = 1.0e-6
u_r = [0.0, 6.28, 48]
v_r = [-0.3, 0.3, 12]
(u, v) = mgrid[u_r[0]:u_r[1]+eps:(u_r[1]-u_r[0])/float(u_r[2]-1),\
                 v_r[0]:v_r[1]+eps:(v_r[1]-v_r[0])/float(v_r[2]-1)]
X = {'ey': v*sin(u)*cos(u/2) + sin(u), 'ex': v*cos(u/2)*cos(u) + cos(u), 'ez': v*sin(u/2)}
scal_tan = 0.15
x = X['ex']
y = X['ey']
z = X['ez']
du = {'ey': -v*sin(u/2)*sin(u)/2 + v*cos(u/2)*cos(u) + cos(u), 'ex': -v*sin(u/2)*cos(u)/2 - v*sin(u)*cos(u/2) - sin(u), 'ez': v*cos(u/2)/2}
dv = {'ey': sin(u)*cos(u/2), 'ex': cos(u/2)*cos(u), 'ez': sin(u/2)}
Zero = zeros(shape(x))
if scal_tan > 0.0:
    du_x = Zero+du['ex']
    du_y = Zero+du['ey']
    du_z = Zero+du['ez']
    dv_x = Zero+dv['ex']
    dv_y = Zero+dv['ey']
    dv_z = Zero+dv['ez']

f = [None]
n = [-v*sin(u/2)**2*sin(u)/2 + v*sin(u/2)*cos(u/2)*cos(u) - v*sin(u)*cos(u/2)**2/2 + sin(u/2)*cos(u), v*sin(u/2)**2*cos(u)/2 + v*sin(u/2)*sin(u)*cos(u/2) + v*cos(u/2)**2*cos(u)/2 + sin(u/2)*sin(u), -v*sin(u)**2*cos(u/2)**2 - v*cos(u/2)**2*cos(u)**2 - sin(u)**2*cos(u/2) - cos(u/2)*cos(u)**2]
skip = [4, 4]
su = skip[0]
sv = skip[1]
if f[0] != None:
    dn_x = f[0]*n[0]
    dn_y = f[0]*n[1]
    dn_z = f[0]*n[2]

from mayavi.mlab import *
figure(bgcolor=(1.0,1.0,1.0))
if False:
    mesh(x,y,z,colormap="gist_earth")
if True:
    for i in range(shape(u)[0]):
        plot3d(x[i,],y[i,],z[i,],line_width=1.0,color=(0.0,0.0,0.0),tube_radius=None)

    xr = swapaxes(x,0,1)
    yr = swapaxes(y,0,1)
    zr = swapaxes(z,0,1)

    for i in range(shape(u)[1]):
        plot3d(xr[i,],yr[i,],zr[i,],line_width=1.0,color=(0.0,0.0,0.0),tube_radius=None)
if scal_tan > 0.0:
    quiver3d(x[::su,::sv],y[::su,::sv],z[::su,::sv],\
             du_x[::su,::sv],du_y[::su,::sv],du_z[::su,::sv],scale_factor=scal_tan,\
             line_width=1.0,color=(0.0,0.0,0.0),scale_mode='vector',mode='arrow',resolution=16)
    quiver3d(x[::su,::sv],y[::su,::sv],z[::su,::sv],\
             dv_x[::su,::sv],dv_y[::su,::sv],dv_z[::su,::sv],scale_factor=scal_tan,\
             line_width=1.0,color=(0.0,0.0,0.0),scale_mode='vector',mode='arrow',resolution=16)
if f[0] != None:
    quiver3d(x[::su,::sv],y[::su,::sv],z[::su,::sv],\
             dn_x[::su,::sv],dn_y[::su,::sv],dn_z[::su,::sv],\
             line_width=1.0,color=(0.0,0.0,0.0),scale_mode='none',mode='cone',\
             resolution=16,opacity=0.5)

