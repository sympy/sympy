from sympy.solvers.deutils import ode_order
from sympy.core.symbol import Dummy


def _get_func_order(eqs, funcs):
    order = dict()
    func_dict = dict()

    for func in funcs:
        if not order.get(func, False):
            max_order = 0
            for i, eq in enumerate(eqs):
                order_ = ode_order(eq, func)
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
    r"""
    Returns a dictionary with details of the eqs if every equation is constant coefficient
    and linear else returns None

    Explanation
    ===========

    This function takes the eqs, converts it into a form Ax = b where x is a vector of terms
    containing dependent variables and their derivatives till their maximum order. If it is
    possible to convert eqs into Ax = b, then all the equations in eqs are linear otherwise
    they are non-linear.

    To check if the equations are constant coefficient, we need to check if all the terms in
    A obtained above are constant or not.

    To check if the equations are homogeneous or not, we need to check if b is a zero matrix
    or not.

    Parameters
    ==========

    eqs: List
        List of ODEs
    funcs: List
        List of dependent variables
    t: Symbol
        Independent variable of the equations in eqs

    Returns
    =======

    match = {
        'no_of_equation': len(eqs),
        'eq': eqs,
        'func': funcs,
        'order': order,
        'is_linear': is_linear,
        'is_constant': is_constant,
        'is_homogeneous': is_homogeneous,
    }

    Dict or None
        Dict with values for keys:
            1. no_of_equation: Number of equations
            2. eq: The set of equations
            3. func: List of dependent variables
            4. order: A dictionary that gives the order of the
                      dependent variable in eqs
            5. is_linear: Boolean value indicating if the set of
                          equations are linear or not.
            6. is_constant: Boolean value indicating if the set of
                          equations have constant coefficients or not.
            7. is_homogeneous: Boolean value indicating if the set of
                          equations are homogeneous or not.
        This Dict is the answer returned if the eqs are linear and constant
        coefficient. Otherwise, None is returned.
    """
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
