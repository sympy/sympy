from sympy.polys.galoispolys import gf_from_dict, gf_factor, gf_factor_sqf

from sympy import pi, nextprime

def gathen_poly(n, p):
    return gf_from_dict({n: 1, 1: 1, 0: 1}, p)

def shoup_poly(n, p):
    f = [1] * (n+1)
    for i in xrange(1, n+1):
        f[i] = (f[i-1]**2 + 1) % p
    return f

def genprime(n):
    return nextprime(int((2**n * pi).evalf()))

p_10 = genprime(10)
f_10 = gathen_poly(10, p_10)

p_20 = genprime(20)
f_20 = gathen_poly(20, p_20)

def timeit_gathen_poly_f10_zassenhaus():
    gf_factor_sqf(f_10, p_10, method='zassenhaus')

def timeit_gathen_poly_f10_shoup():
    gf_factor_sqf(f_10, p_10, method='shoup')

def timeit_gathen_poly_f20_zassenhaus():
    gf_factor_sqf(f_20, p_20, method='zassenhaus')

def timeit_gathen_poly_f20_shoup():
    gf_factor_sqf(f_20, p_20, method='shoup')

P_08  = genprime(8)
F_10 = shoup_poly(10, P_08)

P_18 = genprime(18)
F_20 = shoup_poly(20, P_18)

def timeit_shoup_poly_F10_zassenhaus():
    gf_factor_sqf(F_10, P_08, method='zassenhaus')

def timeit_shoup_poly_F10_shoup():
    gf_factor_sqf(F_10, P_08, method='shoup')

def timeit_shoup_poly_F20_zassenhaus():
    gf_factor_sqf(F_20, P_18, method='zassenhaus')

def timeit_shoup_poly_F20_shoup():
    gf_factor_sqf(F_20, P_18, method='shoup')
