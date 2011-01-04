# Elliptic integral E(z,m) for some different m
f1 = lambda z: ellipe(z,-2)
f2 = lambda z: ellipe(z,-1)
f3 = lambda z: ellipe(z,0)
f4 = lambda z: ellipe(z,1)
f5 = lambda z: ellipe(z,2)
plot([f1,f2,f3,f4,f5], [0,pi], [0,4])
