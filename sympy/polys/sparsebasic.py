"""Basic tools for sparse polynomials in $\K[x]$ or $\K[X]$. """

def smp_validate(f):
    pass

def sup_zero():
    return {}

def smp_zero(u):
    return {}

def smp_one(u, K):
    return {(0,)*(u+1): K.one}

def smp_from_ground(f, u, ord, K):
    pass

def smp_from_dict(f, u, ord, K):
    pass

def smp_from_sympy_dict(f, u, ord, K):
    pass

def smp_from_dict(f, u, ord, K):
    pass

def smp_from_sympy_dict(f, u, ord, K):
    pass

def smp_to_ground(f, u, ord, K):
    pass

def smp_to_dict(f, u, ord, K):
    pass

def smp_to_sympy_dict(f, u, ord, K):
    pass

def smp_to_dict(f, u, ord, K):
    pass

def smp_to_sympy_dict(f, u, ord, K):
    pass

def smp_set_order(f, u, ord, K):
    pass

def smp_set_domain(f, u, ord, K0, K1):
    pass

def smp_ground_LC(f, u, ord, K):
    pass

def smp_ground_LM(f, u, ord, K):
    pass

def smp_ground_LT(f, u, ord, K):
    pass

def smp_ground_TC(f, u, ord, K):
    pass

def smp_ground_TM(f, u, ord, K):
    pass

def smp_ground_TT(f, u, ord, K):
    pass

def smp_ground_EC(f, u, ord, K):
    pass

def smp_ground_EM(f, u, ord, K):
    pass

def smp_ground_ET(f, u, ord, K):
    pass

def smp_coeffs(f, u, ord, K):
    pass

def smp_monoms(f, u, ord, K):
    pass

def smp_terms(f, u, ord, K):
    pass

def smp_all_coeffs(f, u, ord, K):
    pass

def smp_all_monoms(f, u, ord, K):
    pass

def smp_all_terms(f, u, ord, K):
    pass

def smp_degree(f, j, u):
    pass

def smp_degrees(f, u):
    pass

def smp_total_degree(f, u):
    return sum(smp_degrees(f, u))

def smp_deflate(f, u, ord, K):
    pass

def smp_inflate(f, M, u, ord, K):
    pass

def smp_terms_gcd(f, u, ord, K):
    pass
