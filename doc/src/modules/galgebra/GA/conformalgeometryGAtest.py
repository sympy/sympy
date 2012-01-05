    >>> print '\nExample: Conformal representations of circles, lines, spheres, and planes'

    Example: Conformal representations of circles, lines, spheres, and planes

    >>> metric = '1 0 0 0 0,0 1 0 0 0,0 0 1 0 0,0 0 0 0 2,0 0 0 2 0'
    >>> MV.setup('e0 e1 e2 n nbar',metric,debug=0)
    >>> MV.set_str_format(1)
    >>> e = n+nbar

    #conformal representation of points
    >>> A = make_vector(e0)    # point a = (1,0,0)  A = F(a)
    >>> C = make_vector(-1*e0) # point c = (-1,0,0) C = F(c)
    >>> D = make_vector(e2)    # point d = (0,0,1)  D = F(d)
    >>> X = make_vector('x',3)

    >>> print 'a = e0, b = e1, c = -e0, and d = e2'
    a = e0, b = e1, c = -e0, and d = e2

    >>> print 'A = F(a) = 1/2*(a*a*n+2*a-nbar), etc.'
    A = F(a) = 1/2*(a*a*n+2*a-nbar), etc.

    >>> print 'Circle through a, b, and c'
    Circle through a, b, and c

    >>> print 'Circle: A^B^C^X = 0 =',(A^B^C^X)
    Circle: A^B^C^X = 0 = {-x2}e0^e1^e2^n
    +{x2}e0^e1^e2^nbar
    +{-1/2 + 1/2*x0**2 + 1/2*x1**2 + 1/2*x2**2}e0^e1^n^nbar

    >>> print 'Line through a and b'
    Line through a and b
    +{-1/2 + x0/2 + x1/2}e0^e1^n^nbar
    +{x2/2}e0^e2^n^nbar
    +{-x2/2}e1^e2^n^nbar

    >>> print 'Line  : A^B^n^X = 0 =',(A^B^n^X)
    Line  : A^B^n^X = 0 = {-x2}e0^e1^e2^n

    >>> print 'Sphere through a, b, c, and d'
    Sphere through a, b, c, and d

    >>> print 'Sphere: A^B^C^D^X = 0 =',(A^B^C^D^X)
    Sphere: A^B^C^D^X = 0 = {1/2 - 1/2*x0**2 - 1/2*x1**2 - 1/2*x2**2}e0^e1^e2^n^nbar

    >>> print 'Plane through a, b, and d'
    Plane through a, b, and d

    >>> print 'Plane : A^B^n^D^X = 0 =',(A^B^n^D^X)
    Plane : A^B^n^D^X = 0 = {1/2 - x0/2 - x1/2 - x2/2}e0^e1^e2^n^nbar
