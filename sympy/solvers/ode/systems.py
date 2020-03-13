from sympy.solvers.deutils import ode_order
from sympy.core.symbol import Dummy


def _get_func_order(eqs, funcs):
    order = dict()
    func_dict = dict()

    for func in funcs:
        if not order.get(func, False):
            max_order = 0
            for i, eqs_ in enumerate(eqs):
                order_ = ode_order(eqs, func)
                if max_order < order_:
                    max_order = order_
                    eq_no = i
            if eq_no in func_dict:
                list_func = []
                list_func.append(func_dict[eq_no])
                list_func.append(func)
                func_dict[eq_no] = list_func
            else:
                func_dict[eq_no] = func
            order[func] = max_order

    funcs = [func_dict[i] for i in range(len(func_dict))]

    return funcs, order


def neq_nth_linear_constant_coeff_match(eqs, funcs, t):
    from sympy.solvers.solveset import linear_eq_to_matrix

    # Error for i == 0 can be added but isn't for now

    # Assuming funcs has all the funcs mentioned and
    # len(funcs) == len(eqs)

    # Getting the func_dict and order using the helper
    # function
    funcs, order = _get_func_order(eqs, funcs)

    # Not adding the check if the len(func.args) for
    # every func in funcs is 1

    rep = {func.diff(t, n): Dummy() for func in funcs for n in range(order[func] + 1)}
    eqs_sub = [eq.subs(rep) for eq in eqs]

    # Linearity check
    try:
        A, b = linear_eq_to_matrix(eqs_sub, rep.values())

    except ValueError:
        return None

    is_linear = True

    # Constant coefficient check
    is_constant = True
    for coef in A:
        if is_constant == True:
            if coef.as_independent(t, as_Add=True)[1] != 0:
                is_constant = False

    # Homogeneous check
    is_homogeneous = True if b.is_zero_matrix else False

    match = {
        'no_of_equation': len(eqs),
        'eq': eqs,
        'func': funcs,
        'order': order,
        'is_linear': is_linear,
        'is_constant': is_constant,
        'is_homogeneous': is_homogeneous,
    }

    if match['is_linear'] and match['is_constant']:
        return match
    return match
