        print '\nExample: Conformal representations of circles, lines, spheres, and planes'

        metric = '1 0 0 0 0,0 1 0 0 0,0 0 1 0 0,0 0 0 0 2,0 0 0 2 0'

        e0,e1,e2,n,nbar = MV.setup('e0 e1 e2 n nbar',metric,debug=0)
        MV.set_str_format(1)
        e = n+nbar
        #conformal representation of points

        A = F(e0,n,nbar)    # point a = (1,0,0)  A = F(a)
        B = F(e1,n,nbar)    # point b = (0,1,0)  B = F(b)
        C = F(-1*e0,n,nbar) # point c = (-1,0,0) C = F(c)
        D = F(e2,n,nbar)    # point d = (0,0,1)  D = F(d)
        x0,x1,x2 = sympy.symbols('x0 x1 x2')
        X = F(MV([x0,x1,x2],'vector'),n,nbar)

        print 'a = e0, b = e1, c = -e0, and d = e2'
        print 'A = F(a) = 1/2*(a*a*n+2*a-nbar), etc.'
        print 'Circle through a, b, and c'
        print 'Circle: A^B^C^X = 0 =',(A^B^C^X)
        print 'Line through a and b'
        print 'Line  : A^B^n^X = 0 =',(A^B^n^X)
        print 'Sphere through a, b, c, and d'
        print 'Sphere: A^B^C^D^X = 0 =',(A^B^C^D^X)
        print 'Plane through a, b, and d'
        print 'Plane : A^B^n^D^X = 0 =',(A^B^n^D^X)

Example: Conformal representations of circles, lines, spheres, and planes
a = e0, b = e1, c = -e0, and d = e2
A = F(a) = 1/2*(a*a*n+2*a-nbar), etc.
Circle through a, b, and c
Circle: A^B^C^X = 0 = {-x2}e0^e1^e2^n
+{x2}e0^e1^e2^nbar
+{-1/2 + 1/2*x0**2 + 1/2*x1**2 + 1/2*x2**2}e0^e1^n^nbar

Line through a and b
Line  : A^B^n^X = 0 = {-x2}e0^e1^e2^n
+{-1/2 + x0/2 + x1/2}e0^e1^n^nbar
+{x2/2}e0^e2^n^nbar
+{-x2/2}e1^e2^n^nbar

Sphere through a, b, c, and d
Sphere: A^B^C^D^X = 0 = {1/2 - 1/2*x0**2 - 1/2*x1**2 - 1/2*x2**2}e0^e1^e2^n^nbar

Plane through a, b, and d
Plane : A^B^n^D^X = 0 = {1/2 - x0/2 - x1/2 - x2/2}e0^e1^e2^n^nbar
