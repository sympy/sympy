from sympy.external import import_module
matchpy = import_module("matchpy")
if matchpy:
    from matchpy import CustomConstraint
    from sympy.integrals.rubi.utility_function import FreeQ

    constraint_freeq_a = CustomConstraint(lambda a, x: FreeQ(a, x))
    constraint_freeq_b = CustomConstraint(lambda b, x: FreeQ(b, x))
    constraint_freeq_c = CustomConstraint(lambda c, x: FreeQ(c, x))
    constraint_freeq_d = CustomConstraint(lambda d, x: FreeQ(d, x))
    constraint_freeq_e = CustomConstraint(lambda e, x: FreeQ(e, x))
    constraint_freeq_f = CustomConstraint(lambda f, x: FreeQ(f, x))
    constraint_freeq_g = CustomConstraint(lambda g, x: FreeQ(g, x))
    constraint_freeq_h = CustomConstraint(lambda h, x: FreeQ(h, x))
    constraint_freeq_m = CustomConstraint(lambda m, x: FreeQ(m, x))
    constraint_freeq_n = CustomConstraint(lambda n, x: FreeQ(n, x))
    constraint_freeq_p = CustomConstraint(lambda p, x: FreeQ(p, x))
