# Real part of spherical harmonic Y_(4,2)(theta,phi)
def Y(l,m):
    def g(theta,phi):
        R = abs(fp.re(fp.spherharm(l,m,theta,phi)))
        x = R*fp.cos(phi)*fp.sin(theta)
        y = R*fp.sin(phi)*fp.sin(theta)
        z = R*fp.cos(theta)
        return [x,y,z]
    return g

fp.splot(Y(4,2), [0,fp.pi], [0,2*fp.pi], points=300)