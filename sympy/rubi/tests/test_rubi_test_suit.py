from sympy.core.symbol import symbols
from sympy.rubi.rubi import rubi_integrate
from sympy.functions import log, sqrt
from sympy import sympify, N

a, b, c, d, e, f, m, x, u = symbols('a b c d e f m x u')

def test_rubi_algebriac_1_2():
    # Integrands of the form x**m
    assert rubi_integrate(x**m, x) == x**(1 + m)/(1 + m)
    assert rubi_integrate(x**100, x) == 1/101*x**101
    assert rubi_integrate(x**3, x) == 1/4*x**4
    assert rubi_integrate(x**2, x) == 1/3*x**3
    assert rubi_integrate(x, x) == 1/2*x**2
    #assert rubi_integrate(1, x) == x
    assert rubi_integrate(1/x, x) == log(x)
    assert rubi_integrate(1/x**2, x) == (-1)/x
    assert rubi_integrate(1/x**3, x) == (-1/2)/x**2
    assert rubi_integrate(1/x**4, x) == (-1/3)/x**3
    assert rubi_integrate(1/x**100, x) == (-1/99)/x**99

    # Integrands of the form x**(m/2
    assert rubi_integrate(x**(5/2), x) == 2/7*x**(7/2)
    assert rubi_integrate(x**(3/2), x) == 2/5*x**(5/2)
    assert rubi_integrate(x**(1/2), x) == 2/3*x**(3/2)
    assert rubi_integrate(1/x**(1/2), x) == 2*sqrt(x)
    assert rubi_integrate(1/x**(3/2), x) == (-2)/sqrt(x)
    assert rubi_integrate(1/x**(5/2), x) == (-2/3)/x**(3/2)

    # Integrands of the form x**(m/3
    '''
    assert rubi_integrate(x**(5/3), x) == 3/8*x**(8/3)
    assert rubi_integrate(x**(4/3), x) == 3/7*x**(7/3)
    assert rubi_integrate(x**(2/3), x) == 3/5*x**(5/3)
    assert rubi_integrate(x**(1/3), x) == 3/4*x**(4/3)
    assert rubi_integrate(1/x**(1/3), x) == 3/2*x**(2/3)
    assert rubi_integrate(1/x**(2/3), x) == 3*x**(1/3)
    assert rubi_integrate(1/x**(4/3), x) == (-3)/x**(1/3)
    assert rubi_integrate(1/x**(5/3), x) == (-3/2)/x**(2/3)
    '''

    # n>0
    assert rubi_integrate(x**3*(a + b*x), x) == 1/4*a*x**4 + 1/5*b*x**5
    assert rubi_integrate(x**2*(a + b*x), x) == 1/3*a*x**3 + 1/4*b*x**4
    assert rubi_integrate(x*(a + b*x), x) == 1/2*a*x**2 + 1/3*b*x**3
    assert rubi_integrate(a + b*x, x) == a*x + 1/2*b*x**2
    assert rubi_integrate((a + b*x)/x, x) == b*x + a*log(x)
    assert rubi_integrate((a + b*x)/x**2, x) == -a/x + b*log(x)
    assert rubi_integrate((a + b*x)/x**3, x) == -1/2*a/x**2-b/x
    assert rubi_integrate((a + b*x)/x**4, x) == -1/3*a/x**3-1/2*b/x**2
    assert rubi_integrate((a + b*x)/x**5, x) == -1/4*a/x**4-1/3*b/x**3
    assert rubi_integrate(x**3*(a + b*x)**2, x) == 1/4*a**2*x**4 + 2/5*a*b*x**5 + 1/6*b**2*x**6
    assert rubi_integrate(x**2*(a + b*x)**2, x) == 1/3*a**2*x**3 + 1/2*a*b*x**4 + 1/5*b**2*x**5
    assert rubi_integrate(x*(a + b*x)**2, x) == 1/2*a**2*x**2 + 2/3*a*b*x**3 + 1/4*b**2*x**4
    assert rubi_integrate((a + b*x)**2, x) == 1/3*(a + b*x)**3/b
    assert rubi_integrate((a + b*x)**2/x, x) == 2*a*b*x + 1/2*b**2*x**2 + a**2*log(x)
    assert rubi_integrate((a + b*x)**2/x**2, x) == -a**2/x + b**2*x + 2*a*b*log(x)
    assert rubi_integrate((a + b*x)**2/x**3, x) == -1/2*a**2/x**2-2*a*b/x + b**2*log(x)
    assert rubi_integrate((a + b*x)**2/x**4, x) == -1/3*(a + b*x)**3/(a*x**3)
    assert rubi_integrate((a + b*x)**2/x**5, x) == -1/4*a**2/x**4-2/3*a*b/x**3-1/2*b**2/x**2
    assert rubi_integrate((a + b*x)**2/x**6, x) == -1/5*a**2/x**5-1/2*a*b/x**4-1/3*b**2/x**3
    assert rubi_integrate((a + b*x)**2/x**7, x) == -1/6*a**2/x**6-2/5*a*b/x**5-1/4*b**2/x**4
    assert rubi_integrate((a + b*x)**2/x**8, x) == -1/7*a**2/x**7-1/3*a*b/x**6-1/5*b**2/x**5
    assert rubi_integrate(x**3*(a + b*x)**3, x) == 1/4*a**3*x**4 + 3/5*a**2*b*x**5 + 1/2*a*b**2*x**6 + 1/7*b**3*x**7
    assert rubi_integrate(x**2*(a + b*x)**3, x) == 1/3*a**3*x**3 + 3/4*a**2*b*x**4 + 3/5*a*b**2*x**5 + 1/6*b**3*x**6
    assert rubi_integrate(x*(a + b*x)**3, x) == -1/4*a*(a + b*x)**4/b**2 + 1/5*(a + b*x)**5/b**2
    assert rubi_integrate((a + b*x)**3, x) == 1/4*(a + b*x)**4/b
    assert rubi_integrate((a + b*x)**3/x, x) == 3*a**2*b*x + 3/2*a*b**2*x**2 + 1/3*b**3*x**3 + a**3*log(x)
    assert rubi_integrate((a + b*x)**3/x**2, x) == -a**3/x + 3*a*b**2*x + 1/2*b**3*x**2 + 3*a**2*b*log(x)
    assert rubi_integrate((a + b*x)**3/x**3, x) == -1/2*a**3/x**2-3*a**2*b/x + b**3*x + 3*a*b**2*log(x)
    assert rubi_integrate((a + b*x)**3/x**4, x) == -1/3*a**3/x**3-3/2*a**2*b/x**2-3*a*b**2/x + b**3*log(x)
    assert rubi_integrate((a + b*x)**3/x**5, x) == -1/4*(a + b*x)**4/(a*x**4)
    assert rubi_integrate((a + b*x)**3/x**6, x) == -1/5*a**3/x**5-3/4*a**2*b/x**4-a*b**2/x**3-1/2*b**3/x**2
    assert rubi_integrate((a + b*x)**3/x**7, x) == -1/6*a**3/x**6-3/5*a**2*b/x**5-3/4*a*b**2/x**4-1/3*b**3/x**3
    assert rubi_integrate((a + b*x)**3/x**8, x) == -1/7*a**3/x**7-1/2*a**2*b/x**6-3/5*a*b**2/x**5-1/4*b**3/x**4
    assert rubi_integrate(x**7*(a + b*x)**7, x) == 1/8*a**7*x**8 + 7/9*a**6*b*x**9 + 21/10*a**5*b**2*x**10 + 35/11*a**4*b**3*x**11 + 35/12*a**3*b**4*x**12 + 21/13*a**2*b**5*x**13 + 1/2*a*b**6*x**14 + 1/15*b**7*x**15
    assert rubi_integrate(x**6*(a + b*x)**7, x) == 1/7*a**7*x**7 + 7/8*a**6*b*x**8 + 7/3*a**5*b**2*x**9 + 7/2*a**4*b**3*x**10 + 35/11*a**3*b**4*x**11 + 7/4*a**2*b**5*x**12 + 7/13*a*b**6*x**13 + 1/14*b**7*x**14
    assert rubi_integrate(x**5*(a + b*x)**7, x) == 1/6*a**7*x**6 + a**6*b*x**7 + 21/8*a**5*b**2*x**8 + 35/9*a**4*b**3*x**9 + 7/2*a**3*b**4*x**10 + 21/11*a**2*b**5*x**11 + 7/12*a*b**6*x**12 + 1/13*b**7*x**13
    assert rubi_integrate(x**4*(a + b*x)**7, x) == 1/8*a**4*(a + b*x)**8/b**5-4/9*a**3*(a + b*x)**9/b**5 + 3/5*a**2*(a + b*x)**10/b**5-4/11*a*(a + b*x)**11/b**5 + 1/12*(a + b*x)**12/b**5
    assert rubi_integrate(x**3*(a + b*x)**7, x) == -1/8*a**3*(a + b*x)**8/b**4 + 1/3*a**2*(a + b*x)**9/b**4-3/10*a*(a + b*x)**10/b**4 + 1/11*(a + b*x)**11/b**4
    assert rubi_integrate(x**2*(a + b*x)**7, x) == 1/8*a**2*(a + b*x)**8/b**3-2/9*a*(a + b*x)**9/b**3 + 1/10*(a + b*x)**10/b**3
    assert rubi_integrate(x*(a + b*x)**7, x) == -1/8*a*(a + b*x)**8/b**2 + 1/9*(a + b*x)**9/b**2
    assert rubi_integrate((a + b*x)**7, x) == 1/8*(a + b*x)**8/b
    assert rubi_integrate((a + b*x)**7/x, x) == 7*a**6*b*x + 21/2*a**5*b**2*x**2 + 35/3*a**4*b**3*x**3 + 35/4*a**3*b**4*x**4 + 21/5*a**2*b**5*x**5 + 7/6*a*b**6*x**6 + 1/7*b**7*x**7 + a**7*log(x)
    assert rubi_integrate((a + b*x)**7/x**2, x) == -a**7/x + 21*a**5*b**2*x + 35/2*a**4*b**3*x**2 + 35/3*a**3*b**4*x**3 + 21/4*a**2*b**5*x**4 + 7/5*a*b**6*x**5 + 1/6*b**7*x**6 + 7*a**6*b*log(x)
    assert rubi_integrate((a + b*x)**7/x**3, x) == -1/2*a**7/x**2-7*a**6*b/x + 35*a**4*b**3*x + 35/2*a**3*b**4*x**2 + 7*a**2*b**5*x**3 + 7/4*a*b**6*x**4 + 1/5*b**7*x**5 + 21*a**5*b**2*log(x)
    assert rubi_integrate((a + b*x)**7/x**4, x) == -1/3*a**7/x**3-7/2*a**6*b/x**2-21*a**5*b**2/x + 35*a**3*b**4*x + 21/2*a**2*b**5*x**2 + 7/3*a*b**6*x**3 + 1/4*b**7*x**4 + 35*a**4*b**3*log(x)
    assert rubi_integrate((a + b*x)**7/x**5, x) == -1/4*a**7/x**4-7/3*a**6*b/x**3-21/2*a**5*b**2/x**2-35*a**4*b**3/x + 21*a**2*b**5*x + 7/2*a*b**6*x**2 + 1/3*b**7*x**3 + 35*a**3*b**4*log(x)
    assert rubi_integrate((a + b*x)**7/x**6, x) == -1/5*a**7/x**5-7/4*a**6*b/x**4-7*a**5*b**2/x**3-35/2*a**4*b**3/x**2-35*a**3*b**4/x + 7*a*b**6*x + 1/2*b**7*x**2 + 21*a**2*b**5*log(x)
    assert rubi_integrate((a + b*x)**7/x**7, x) == -1/6*a**7/x**6-7/5*a**6*b/x**5-21/4*a**5*b**2/x**4-35/3*a**4*b**3/x**3-35/2*a**3*b**4/x**2-21*a**2*b**5/x + b**7*x + 7*a*b**6*log(x)
    assert rubi_integrate((a + b*x)**7/x**8, x) == -1/7*a**7/x**7-7/6*a**6*b/x**6-21/5*a**5*b**2/x**5-35/4*a**4*b**3/x**4-35/3*a**3*b**4/x**3-21/2*a**2*b**5/x**2-7*a*b**6/x + b**7*log(x)
    assert rubi_integrate((a + b*x)**7/x**9, x) == -1/8*(a + b*x)**8/(a*x**8)
    assert rubi_integrate((a + b*x)**7/x**10, x) == -1/9*(a + b*x)**8/(a*x**9) + 1/72*b*(a + b*x)**8/(a**2*x**8)
    assert rubi_integrate((a + b*x)**7/x**11, x) == -1/10*(a + b*x)**8/(a*x**10) + 1/45*b*(a + b*x)**8/(a**2*x**9)-1/360*b**2*(a + b*x)**8/(a**3*x**8)
    assert rubi_integrate((a + b*x)**7/x**12, x) == -1/11*(a + b*x)**8/(a*x**11) + 3/110*b*(a + b*x)**8/(a**2*x**10)-1/165*b**2*(a + b*x)**8/(a**3*x**9) + 1/1320*b**3*(a + b*x)**8/(a**4*x**8)
    assert rubi_integrate((a + b*x)**7/x**13, x) == -1/12*a**7/x**12-7/11*a**6*b/x**11-21/10*a**5*b**2/x**10-35/9*a**4*b**3/x**9-35/8*a**3*b**4/x**8-3*a**2*b**5/x**7-7/6*a*b**6/x**6-1/5*b**7/x**5
    assert rubi_integrate((a + b*x)**7/x**14, x) == -1/13*a**7/x**13-7/12*a**6*b/x**12-21/11*a**5*b**2/x**11-7/2*a**4*b**3/x**10-35/9*a**3*b**4/x**9-21/8*a**2*b**5/x**8-a*b**6/x**7-1/6*b**7/x**6
    assert rubi_integrate((a + b*x)**7/x**15, x) == -1/14*a**7/x**14-7/13*a**6*b/x**13-7/4*a**5*b**2/x**12-35/11*a**4*b**3/x**11-7/2*a**3*b**4/x**10-7/3*a**2*b**5/x**9-7/8*a*b**6/x**8-1/7*b**7/x**7
    assert rubi_integrate((a + b*x)**7/x**16, x) == -1/15*a**7/x**15-1/2*a**6*b/x**14-21/13*a**5*b**2/x**13-35/12*a**4*b**3/x**12-35/11*a**3*b**4/x**11-21/10*a**2*b**5/x**10-7/9*a*b**6/x**9-1/8*b**7/x**8
    assert rubi_integrate(x**9*(a + b*x)**10, x) == 1/10*a**10*x**10 + 10/11*a**9*b*x**11 + 15/4*a**8*b**2*x**12 + 120/13*a**7*b**3*x**13 + 15*a**6*b**4*x**14 + 84/5*a**5*b**5*x**15 + 105/8*a**4*b**6*x**16 + 120/17*a**3*b**7*x**17 + 5/2*a**2*b**8*x**18 + 10/19*a*b**9*x**19 + 1/20*b**10*x**20
    assert rubi_integrate(x**8*(a + b*x)**10, x) == 1/9*a**10*x**9 + a**9*b*x**10 + 45/11*a**8*b**2*x**11 + 10*a**7*b**3*x**12 + 210/13*a**6*b**4*x**13 + 18*a**5*b**5*x**14 + 14*a**4*b**6*x**15 + 15/2*a**3*b**7*x**16 + 45/17*a**2*b**8*x**17 + 5/9*a*b**9*x**18 + 1/19*b**10*x**19
    assert rubi_integrate(x**7*(a + b*x)**10, x) == 1/8*a**10*x**8 + 10/9*a**9*b*x**9 + 9/2*a**8*b**2*x**10 + 120/11*a**7*b**3*x**11 + 35/2*a**6*b**4*x**12 + 252/13*a**5*b**5*x**13 + 15*a**4*b**6*x**14 + 8*a**3*b**7*x**15 + 45/16*a**2*b**8*x**16 + 10/17*a*b**9*x**17 + 1/18*b**10*x**18
    assert rubi_integrate(x**6*(a + b*x)**10, x) == 1/11*a**6*(a + b*x)**11/b**7-1/2*a**5*(a + b*x)**12/b**7 + 15/13*a**4*(a + b*x)**13/b**7-10/7*a**3*(a + b*x)**14/b**7 + a**2*(a + b*x)**15/b**7-3/8*a*(a + b*x)**16/b**7 + 1/17*(a + b*x)**17/b**7
    assert rubi_integrate(x**5*(a + b*x)**10, x) == -1/11*a**5*(a + b*x)**11/b**6 + 5/12*a**4*(a + b*x)**12/b**6-10/13*a**3*(a + b*x)**13/b**6 + 5/7*a**2*(a + b*x)**14/b**6-1/3*a*(a + b*x)**15/b**6 + 1/16*(a + b*x)**16/b**6
    assert rubi_integrate(x**4*(a + b*x)**10, x) == 1/11*a**4*(a + b*x)**11/b**5-1/3*a**3*(a + b*x)**12/b**5 + 6/13*a**2*(a + b*x)**13/b**5-2/7*a*(a + b*x)**14/b**5 + 1/15*(a + b*x)**15/b**5
    assert rubi_integrate(x**3*(a + b*x)**10, x) == -1/11*a**3*(a + b*x)**11/b**4 + 1/4*a**2*(a + b*x)**12/b**4-3/13*a*(a + b*x)**13/b**4 + 1/14*(a + b*x)**14/b**4
    assert rubi_integrate(x**2*(a + b*x)**10, x) == 1/11*a**2*(a + b*x)**11/b**3-1/6*a*(a + b*x)**12/b**3 + 1/13*(a + b*x)**13/b**3
    assert rubi_integrate(x*(a + b*x)**10, x) == -1/11*a*(a + b*x)**11/b**2 + 1/12*(a + b*x)**12/b**2
    assert rubi_integrate((a + b*x)**10, x) == 1/11*(a + b*x)**11/b
    assert rubi_integrate((a + b*x)**10/x, x) == 10*a**9*b*x + 45/2*a**8*b**2*x**2 + 40*a**7*b**3*x**3 + 105/2*a**6*b**4*x**4 + 252/5*a**5*b**5*x**5 + 35*a**4*b**6*x**6 + 120/7*a**3*b**7*x**7 + 45/8*a**2*b**8*x**8 + 10/9*a*b**9*x**9 + 1/10*b**10*x**10 + a**10*log(x)
    assert rubi_integrate((a + b*x)**10/x**2, x) == -a**10/x + 45*a**8*b**2*x + 60*a**7*b**3*x**2 + 70*a**6*b**4*x**3 + 63*a**5*b**5*x**4 + 42*a**4*b**6*x**5 + 20*a**3*b**7*x**6 + 45/7*a**2*b**8*x**7 + 5/4*a*b**9*x**8 + 1/9*b**10*x**9 + 10*a**9*b*log(x)
    assert rubi_integrate((a + b*x)**10/x**3, x) == -1/2*a**10/x**2-10*a**9*b/x + 120*a**7*b**3*x + 105*a**6*b**4*x**2 + 84*a**5*b**5*x**3 + 105/2*a**4*b**6*x**4 + 24*a**3*b**7*x**5 + 15/2*a**2*b**8*x**6 + 10/7*a*b**9*x**7 + 1/8*b**10*x**8 + 45*a**8*b**2*log(x)
    assert rubi_integrate((a + b*x)**10/x**4, x) == -1/3*a**10/x**3-5*a**9*b/x**2-45*a**8*b**2/x + 210*a**6*b**4*x + 126*a**5*b**5*x**2 + 70*a**4*b**6*x**3 + 30*a**3*b**7*x**4 + 9*a**2*b**8*x**5 + 5/3*a*b**9*x**6 + 1/7*b**10*x**7 + 120*a**7*b**3*log(x)
    assert rubi_integrate((a + b*x)**10/x**5, x) == -1/4*a**10/x**4-10/3*a**9*b/x**3-45/2*a**8*b**2/x**2-120*a**7*b**3/x + 252*a**5*b**5*x + 105*a**4*b**6*x**2 + 40*a**3*b**7*x**3 + 45/4*a**2*b**8*x**4 + 2*a*b**9*x**5 + 1/6*b**10*x**6 + 210*a**6*b**4*log(x)
    assert rubi_integrate((a + b*x)**10/x**6, x) == -1/5*a**10/x**5-5/2*a**9*b/x**4-15*a**8*b**2/x**3-60*a**7*b**3/x**2-210*a**6*b**4/x + 210*a**4*b**6*x + 60*a**3*b**7*x**2 + 15*a**2*b**8*x**3 + 5/2*a*b**9*x**4 + 1/5*b**10*x**5 + 252*a**5*b**5*log(x)
    assert rubi_integrate((a + b*x)**10/x**7, x) == -1/6*a**10/x**6-2*a**9*b/x**5-45/4*a**8*b**2/x**4-40*a**7*b**3/x**3-105*a**6*b**4/x**2-252*a**5*b**5/x + 120*a**3*b**7*x + 45/2*a**2*b**8*x**2 + 10/3*a*b**9*x**3 + 1/4*b**10*x**4 + 210*a**4*b**6*log(x)
    assert rubi_integrate((a + b*x)**10/x**8, x) == -1/7*a**10/x**7-5/3*a**9*b/x**6-9*a**8*b**2/x**5-30*a**7*b**3/x**4-70*a**6*b**4/x**3-126*a**5*b**5/x**2-210*a**4*b**6/x + 45*a**2*b**8*x + 5*a*b**9*x**2 + 1/3*b**10*x**3 + 120*a**3*b**7*log(x)
    assert rubi_integrate((a + b*x)**10/x**9, x) == -1/8*a**10/x**8-10/7*a**9*b/x**7-15/2*a**8*b**2/x**6-24*a**7*b**3/x**5-105/2*a**6*b**4/x**4-84*a**5*b**5/x**3-105*a**4*b**6/x**2-120*a**3*b**7/x + 10*a*b**9*x + 1/2*b**10*x**2 + 45*a**2*b**8*log(x)
    assert rubi_integrate((a + b*x)**10/x**10, x) == -1/9*a**10/x**9-5/4*a**9*b/x**8-45/7*a**8*b**2/x**7-20*a**7*b**3/x**6-42*a**6*b**4/x**5-63*a**5*b**5/x**4-70*a**4*b**6/x**3-60*a**3*b**7/x**2-45*a**2*b**8/x + b**10*x + 10*a*b**9*log(x)
    assert rubi_integrate((a + b*x)**10/x**11, x) == -1/10*a**10/x**10-10/9*a**9*b/x**9-45/8*a**8*b**2/x**8-120/7*a**7*b**3/x**7-35*a**6*b**4/x**6-252/5*a**5*b**5/x**5-105/2*a**4*b**6/x**4-40*a**3*b**7/x**3-45/2*a**2*b**8/x**2-10*a*b**9/x + b**10*log(x)
    assert rubi_integrate((a + b*x)**10/x**12, x) == -1/11*(a + b*x)**11/(a*x**11)
    assert rubi_integrate((a + b*x)**10/x**13, x) == -1/12*(a + b*x)**11/(a*x**12) + 1/132*b*(a + b*x)**11/(a**2*x**11)
    assert rubi_integrate((a + b*x)**10/x**14, x) == -1/13*(a + b*x)**11/(a*x**13) + 1/78*b*(a + b*x)**11/(a**2*x**12)-1/858*b**2*(a + b*x)**11/(a**3*x**11)
    assert rubi_integrate((a + b*x)**10/x**15, x) == -1/14*(a + b*x)**11/(a*x**14) + 3/182*b*(a + b*x)**11/(a**2*x**13)-1/364*b**2*(a + b*x)**11/(a**3*x**12) + 1/4004*b**3*(a + b*x)**11/(a**4*x**11)
    assert rubi_integrate((a + b*x)**10/x**16, x) == -1/15*(a + b*x)**11/(a*x**15) + 2/105*b*(a + b*x)**11/(a**2*x**14)-2/455*b**2*(a + b*x)**11/(a**3*x**13) + 1/1365*b**3*(a + b*x)**11/(a**4*x**12)-1/15015*b**4*(a + b*x)**11/(a**5*x**11)
    assert rubi_integrate((a + b*x)**10/x**17, x) == -1/16*(a + b*x)**11/(a*x**16) + 1/48*b*(a + b*x)**11/(a**2*x**15)-1/168*b**2*(a + b*x)**11/(a**3*x**14) + 1/728*b**3*(a + b*x)**11/(a**4*x**13)-1/4368*b**4*(a + b*x)**11/(a**5*x**12) + 1/48048*b**5*(a + b*x)**11/(a**6*x**11)
    assert rubi_integrate((a + b*x)**10/x**18, x) == -1/17*a**10/x**17-5/8*a**9*b/x**16-3*a**8*b**2/x**15-60/7*a**7*b**3/x**14-210/13*a**6*b**4/x**13-21*a**5*b**5/x**12-210/11*a**4*b**6/x**11-12*a**3*b**7/x**10-5*a**2*b**8/x**9-5/4*a*b**9/x**8-1/7*b**10/x**7
    assert rubi_integrate((a + b*x)**10/x**19, x) == -1/18*a**10/x**18-10/17*a**9*b/x**17-45/16*a**8*b**2/x**16-8*a**7*b**3/x**15-15*a**6*b**4/x**14-252/13*a**5*b**5/x**13-35/2*a**4*b**6/x**12-120/11*a**3*b**7/x**11-9/2*a**2*b**8/x**10-10/9*a*b**9/x**9-1/8*b**10/x**8
    assert rubi_integrate((a + b*x)**10/x**20, x) == -1/19*a**10/x**19-5/9*a**9*b/x**18-45/17*a**8*b**2/x**17-15/2*a**7*b**3/x**16-14*a**6*b**4/x**15-18*a**5*b**5/x**14-210/13*a**4*b**6/x**13-10*a**3*b**7/x**12-45/11*a**2*b**8/x**11-a*b**9/x**10-1/9*b**10/x**9
    assert rubi_integrate((c + d)*(a + b*x)/e, x) == 1/2*(c + d)*(a + b*x)**2/(b*e)
