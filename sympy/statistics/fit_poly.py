from sympy import Float, Matrix, S, Symbol, Integer, lambdify, sqrt

def fit_poly(x, y, sym, degree=1):

    """
    Fit polynomial of degree k through n points (x[i],y[i]), i=0..n
    Returns tuple (fit function, residual)
    Implementation: Least squares of vertical offset.

    Example: Fit line
    ---------------------
    x = Symbol('x')
    p,r = fit_poly([0,8], [1,5], x, 1)
    p
    >>> x/2 + 1
    r
    >>> 0

    Example: Fit parabula
    ---------------------
    data_x = [0,1,2,5]
    data_y = [1,1,3,5]
    x = Symbol('x')
    p,r = fit_poly(data_x, data_y, x, 2)
    p
    >>> -3/181*x**2 + 171/181*x + 133/181
    r
    >>> 0.840941329964219
    """

    def fit_matrix(x,y, degree):
        """calculate matrix for fit"""
        n = len(x)
        k = degree+1
        array=[ [0 for i in range(k)] for j in range(n)]
        for j in range(n):
            for i in range(k):
                array[j][i] = x[j]**i
        return Matrix(array)

    def fit_coeffs(m, v):
        """calculate polynomial coefficients"""
        mt = m.transpose()
        return list( (mt*m).inv()*mt*Matrix(v) )
    
    def residual(x, y, e):
        """calculate residual"""
        f = lambdify( e.atoms(Symbol).pop(), e )
        r = Float(0.0)
        for i in range(len(x)):
            r += ( y[i] - f(x[i]) )**2
        return sqrt(r)
    
    def expr_from_coeff(coeff, sym):
        e = coeff[0]
        for i in range(1,len(coeff)):
            e += coeff[i]*sym**i
        return e

    # check and normalize arguments
    assert( isinstance(x, (tuple,list)) )
    assert( isinstance(y, (tuple,list)) )
    assert( isinstance(sym, Symbol) )
    assert( len(x) == len(y) )
    degree = int(degree)
    assert( int(0) <= degree )
    assert( len(x) > degree )
    for i in range(len(x)):
        x[i] = S(x[i])
        y[i] = S(y[i])
        if not x[i].is_Number:
            x[i] = x[i].n()
        if not y[i].is_Number:
            y[i] = y[i].n()

    # fit polynomial and calculate residual
    m = fit_matrix(x, y, degree)
    coeff = fit_coeffs(m, y)
    e = expr_from_coeff(coeff, sym)
    res = Integer(0)
    if degree != len(x)-1:
        res = residual(x, y, e)
    return e, res

