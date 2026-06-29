from __future__ import annotations

from sympy import (
    Expr, Symbol, cancel, solve, S, factorial, harmonic, degree
)
from sympy.concrete.karr_field import KarrField, KarrElement

def build_karr_field(f: Expr, var: Symbol):
    """Automatically constructs a Karr difference field tower for a given summand."""
    field = KarrField()
    field.add_base(var)

    found_pis = []
    found_sigmas = []

    def traverse(expr):
        if expr == var:
            return
        if expr.is_Pow and expr.exp == var:
            found_pis.append(('pow', expr.base, expr))
        elif isinstance(expr, factorial) and expr.args[0] == var:
            found_pis.append(('factorial', None, expr))
        elif isinstance(expr, harmonic) and expr.args[0] == var:
            found_sigmas.append(('harmonic', None, expr))
        else:
            for arg in expr.args:
                traverse(arg)

    traverse(f)

    expr_map = {}
    gen_idx = 1

    for pi_type, param, expr in found_pis:
        if expr in expr_map:
            continue
        t = Symbol(f"t{gen_idx}")
        gen_idx += 1
        expr_map[expr] = t
        if pi_type == 'pow':
            field.add_pi(t, param)
        elif pi_type == 'factorial':
            field.add_pi(t, var + S.One)

    for sig_type, param, expr in found_sigmas:
        if expr in expr_map:
            continue
        t = Symbol(f"t{gen_idx}")
        gen_idx += 1
        expr_map[expr] = t
        if sig_type == 'harmonic':
            field.add_sigma(t, S.One / (var + S.One))

    f_in_field = f.subs(expr_map)
    return field, f_in_field, expr_map


def solve_difference_equation(a: Expr, f: Expr, field: KarrField) -> Expr:
    """Solves the difference equation sigma(g) - a * g = f over the difference field."""
    num_f, den_f = f.as_numer_denom()

    def solve_undetermined(Q: Expr) -> Expr:
        # 1. Determine degree bounds for each generator in the tower
        deg_bounds = {}
        for t in field.ordered_generators:
            deg_f = degree(num_f, t)
            if deg_f < 0:
                deg_f = 0

            g_type, param = field.generators[t]
            if g_type == 'base' or g_type == 'sigma':
                if a == S.One:
                    deg_bounds[t] = deg_f + 1
                else:
                    deg_bounds[t] = deg_f
            elif g_type == 'pi':
                deg_bounds[t] = deg_f + 1
            else:
                deg_bounds[t] = deg_f

        # 2. Construct the general polynomial P with undetermined coefficients
        monomials = [S.One]
        for t in field.ordered_generators:
            new_monomials = []
            bound = deg_bounds[t]
            for m in monomials:
                for deg in range(bound + 1):
                    new_monomials.append(m * t**deg)
            monomials = new_monomials

        coeff_symbols = []
        P_expr = S.Zero
        for i, m in enumerate(monomials):
            c = Symbol(f"c{i}")
            coeff_symbols.append(c)
            P_expr += c * m

        # 3. Compute the residual = sigma(P/Q) - a * (P/Q) - f
        g_elem = field.element(P_expr / Q)
        sigma_g = g_elem.shift().expr
        res = cancel(sigma_g - a * (P_expr / Q) - f)

        # 4. Clear denominator and solve for undetermined coefficients
        num, den = res.as_numer_denom()

        eqs = []
        def extract_coeffs(expr, vars_list):
            if not vars_list:
                if expr != S.Zero:
                    eqs.append(expr)
                return
            v = vars_list[0]
            deg = degree(expr, v)
            for d in range(deg + 1):
                coeff = expr.coeff(v, d)
                coeff = coeff.subs(v, S.Zero)
                extract_coeffs(coeff, vars_list[1:])

        extract_coeffs(num, field.ordered_generators)

        sol = solve(eqs, coeff_symbols)
        if not sol:
            raise ValueError("No solution exists for this denominator.")

        final_P = P_expr.subs(sol)
        final_P = final_P.subs({c: S.Zero for c in coeff_symbols})
        return cancel(final_P / Q)

    # First, try to solve assuming g is a polynomial (Q = 1)
    from sympy.polys.polyerrors import PolynomialError
    try:
        return solve_undetermined(S.One)
    except (ValueError, PolynomialError):
        pass

    # If that fails, and f has a denominator, try to construct a rational solution
    if den_f == S.One:
        raise ValueError("No solution exists.")

    from sympy import factor
    factored_den = factor(den_f)
    factors = []
    if factored_den.is_Mul:
        factors = list(factored_den.args)
    else:
        factors = [factored_den]

    # Try candidate denominators Q
    for q_candidate in factors:
        if q_candidate.is_number:
            continue
        try:
            return solve_undetermined(q_candidate)
        except (ValueError, PolynomialError):
            pass

    raise ValueError("No rational solution exists.")


def karr_sum(f: Expr, limit_tuple: tuple) -> Expr:
    """Computes the symbolic sum of f over the limit_tuple (var, start, end) using Karr's algorithm."""
    var, start, end = limit_tuple
    field, f_in_field, expr_map = build_karr_field(f, var)

    # Solve sigma(g) - g = f_in_field (a = 1)
    try:
        g_in_field = solve_difference_equation(S.One, f_in_field, field)
    except ValueError:
        raise ValueError(f"Summation of {f} could not be solved by Karr's algorithm.")

    # Map the solution back to the original SymPy expressions
    reverse_map = {v: k for k, v in expr_map.items()}
    g_expr = cancel(g_in_field.subs(reverse_map))

    # The sum is g(end + 1) - g(start)
    sum_expr = cancel(g_expr.subs(var, end + S.One) - g_expr.subs(var, start))
    return sum_expr
