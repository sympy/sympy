from sympy import symbols, Poly, solve_undetermined_coeffs, series
from sympy.abc import x

def pade(f, M, N):
    
    '''
    Calculates the pade approximant p/q to the function f at the origin, where deg(p) = M and deg(q) = N
    p = p_0 + ... + p_Mx^M  and  q = 1 + q_1x + ... + q_Nx^N
    
    
    Example
    ========
    
    >>> from sympy.series import pade
    >>> from sympy import exp
    >>> pade(exp(x), 2, 2)
        {p_0: 1, p_1: 1/2, p_2: 1/12, q_1: -1/2, q_2: 1/12}         
    >>>pade(exp(x), 3, 3)
       {p_0: 1, p_1: 1/2, p_2: 1/10, p_3: 1/120, q_1: -1/2, q_2: 1/10, q_3: -1/120}
    '''
    
    
    '''
    The coefficients are calculated and returned to the user
    
    I want to extract the coefficients from solve_undetermined_coeffs(t, coeffs) and then create the polynomials p and q and then return (p,q) to the user, but I don't know how to extract the p_i and the q_j from solve_undetermined_coeffs(t, coeffs)
    I saw on the net that if y = solve_undetermined_coeffs() and y = {a: 1, b: 1/2, c: 1/12} then the values can be accessed by doing c_0 = y[a]
    But I tried v_2 = u[p_2] and it didn't work
    '''    
        
    p = symbols(f"p_0:{M+1}")

    p_dict = {f"p_{i}": p[i] for i in range (M+1)}
                
    locals().update(p_dict)

    q = symbols(f"q_0:{N+1}")

    q_dict = {f"q_{i}": q[i] for i in range (N+1)}
                
    locals().update(q_dict)

    e = [1]

    for i in range (1, N + 1):
        e.append(q[i])


    p = list(reversed(p))

    q = list(reversed(e))

    
    r = Poly.from_list(p,x)

    r = r.as_expr()

    s = Poly.from_list(q,x)

    s = s.as_expr()

    g = f*s - r

    h = series(g, x, 0, M + N + 1)

    t = h.removeO()

    q = q[:-1]

    coeffs = p + q

    u = solve_undetermined_coeffs(t, coeffs)

    return u  
    

    



    

            
