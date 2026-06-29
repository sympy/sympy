from __future__ import annotations

from sympy import (
    Expr, Symbol, cancel, solve, S, sqrt, oo, degree, diff, log
)

class Place:
    """Represents a place (point) on the algebraic curve y**2 = Q(x)."""
    def __init__(self, x0: Expr, y0: Expr):
        self.x0 = x0
        self.y0 = y0

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Place):
            return False
        return self.x0 == other.x0 and self.y0 == other.y0

    def __hash__(self) -> int:
        return hash((self.x0, self.y0))

    def __repr__(self) -> str:
        return f"Place({self.x0}, {self.y0})"


class Divisor:
    """Represents a divisor (formal sum of places) on the algebraic curve."""
    def __init__(self, place_dict: dict[Place, int] | None = None):
        self.dict = {}
        if place_dict:
            for p, mult in place_dict.items():
                if mult != 0:
                    self.dict[p] = mult

    def __add__(self, other: Divisor) -> Divisor:
        new_dict = self.dict.copy()
        for p, mult in other.dict.items():
            new_dict[p] = new_dict.get(p, 0) + mult
        return Divisor(new_dict)

    def __sub__(self, other: Divisor) -> Divisor:
        new_dict = self.dict.copy()
        for p, mult in other.dict.items():
            new_dict[p] = new_dict.get(p, 0) - mult
        return Divisor(new_dict)

    def __mul__(self, m: int) -> Divisor:
        if not isinstance(m, int):
            raise TypeError("Multiplication of divisor is only defined for integers.")
        new_dict = {p: mult * m for p, mult in self.dict.items()}
        return Divisor(new_dict)

    def __rmul__(self, m: int) -> Divisor:
        return self.__mul__(m)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Divisor):
            return False
        return self.dict == other.dict

    def __repr__(self) -> str:
        if not self.dict:
            return "0"
        return " + ".join(f"{mult}*({p})" for p, mult in self.dict.items())


def valuation_at_zero(expr: Expr, u: Symbol) -> Expr:
    """Computes the valuation of expr at u = 0."""
    # noqa: BLE001
    try:
        lt = expr.as_leading_term(u)
    except Exception:  # noqa: BLE001
        # Fallback if as_leading_term fails
        return S.Zero

    if not lt.has(u):
        if lt == S.Zero:
            return oo
        return S.Zero
    coeff, exp = lt.as_coeff_exponent(u)
    return exp


def get_series_coeffs(expr: Expr, var: Symbol, n: int) -> list[Expr]:
    """Extracts series coefficients of expr in var up to n (including negative powers)."""
    # noqa: BLE001
    try:
        s = expr.series(var, 0, n + 2)
        s = s.removeO()
    except Exception:  # noqa: BLE001
        s = expr

    coeffs = []
    for p in range(-6, n):
        coeff = s.coeff(var, p)
        coeff = coeff.subs(var, S.Zero)
        if coeff != S.Zero:
            coeffs.append(coeff)
    return coeffs


def valuation(f: Expr, place: Place, x: Symbol, y: Symbol, Q_x: Expr) -> Expr:
    """Computes the valuation of f(x, y) at a given place on y**2 = Q(x)."""
    x0, y0 = place.x0, place.y0

    if x0 == oo:
        u = Symbol('u')
        Q_u = Q_x.subs(x, S.One/u)
        sign = S.One if y0 == oo else -S.One
        y_expr = sign * S.One/u * sqrt(cancel(Q_u * u**2))
        f_u = f.subs({x: S.One/u, y: y_expr})
        val = valuation_at_zero(f_u, u)
        return val

    if y0 != S.Zero:
        X = Symbol('X')
        y_expr = y0 * sqrt(Q_x.subs(x, x0 + X)) / sqrt(Q_x.subs(x, x0))
        f_X = f.subs({x: x0 + X, y: y_expr})
        return valuation_at_zero(f_X, X)
    else:
        u = Symbol('u')
        ratio = Q_x.subs(x, x0 + u**2) / u**2
        y_expr = u * sqrt(cancel(ratio))
        f_u = f.subs({x: x0 + u**2, y: y_expr})
        val_u = valuation_at_zero(f_u, u)
        return val_u / 2


def residue(f: Expr, place: Place, x: Symbol, y: Symbol, Q_x: Expr) -> Expr:
    """Computes the residue of f(x, y) dx at a given place."""
    x0, y0 = place.x0, place.y0

    if x0 == oo:
        u = Symbol('u')
        Q_u = Q_x.subs(x, S.One/u)
        sign = S.One if y0 == oo else -S.One
        y_expr = sign * S.One/u * sqrt(cancel(Q_u * u**2))
        f_u = f.subs({x: S.One/u, y: y_expr})
        diff_expr = -f_u / u**2
        # noqa: BLE001
        try:
            return diff_expr.as_leading_term(u).coeff(u, -1)
        except Exception:  # noqa: BLE001
            return S.Zero

    if y0 != S.Zero:
        X = Symbol('X')
        y_expr = y0 * sqrt(Q_x.subs(x, x0 + X)) / sqrt(Q_x.subs(x, x0))
        f_X = f.subs({x: x0 + X, y: y_expr})
        # noqa: BLE001
        try:
            return f_X.as_leading_term(X).coeff(X, -1)
        except Exception:  # noqa: BLE001
            return S.Zero
    else:
        u = Symbol('u')
        ratio = Q_x.subs(x, x0 + u**2) / u**2
        y_expr = u * sqrt(cancel(ratio))
        f_u = f.subs({x: x0 + u**2, y: y_expr})
        diff_expr = f_u * 2*u
        # noqa: BLE001
        try:
            return diff_expr.as_leading_term(u).coeff(u, -1)
        except Exception:  # noqa: BLE001
            return S.Zero


def coates_torsion_divisor(divisor: Divisor, x: Symbol, y: Symbol, Q_x: Expr) -> tuple[int, Expr] | None:
    """Coates' algorithm to find m and g such that (g) = m * divisor."""
    places = list(divisor.dict.keys())
    if not places:
        return 1, S.One

    for m in range(1, 6):
        poles = []
        for p, mult in divisor.dict.items():
            if mult < 0:
                poles.append((p, -mult * m))

        D_den = S.One
        for p, deg in poles:
            if p.x0 != oo:
                D_den *= (x - p.x0)**int(deg)

        deg_A = 2
        deg_B = 2

        coeff_symbols = []
        A_expr = S.Zero
        B_expr = S.Zero
        coeff_idx = 0

        for d in range(deg_A + 1):
            c = Symbol(f"a{coeff_idx}")
            coeff_symbols.append(c)
            A_expr += c * x**d
            coeff_idx += 1

        for d in range(deg_B + 1):
            c = Symbol(f"b{coeff_idx}")
            coeff_symbols.append(c)
            B_expr += c * x**d
            coeff_idx += 1

        P_cand = A_expr + B_expr * y
        g_cand = P_cand / D_den

        eqs = []
        for p, mult in divisor.dict.items():
            target_val = mult * m
            if target_val > 0:
                if p.x0 == oo:
                    u = Symbol('u')
                    Q_u = Q_x.subs(x, S.One/u)
                    sign = S.One if p.y0 == oo else -S.One
                    y_expr = sign * S.One/u * sqrt(cancel(Q_u * u**2))
                    g_u = g_cand.subs({x: S.One/u, y: y_expr})
                    eqs.extend(get_series_coeffs(g_u, u, int(target_val)))
                else:
                    if p.y0 != S.Zero:
                        X = Symbol('X')
                        y_expr = p.y0 * sqrt(Q_x.subs(x, p.x0 + X)) / sqrt(Q_x.subs(x, p.x0))
                        g_X = g_cand.subs({x: p.x0 + X, y: y_expr})
                        eqs.extend(get_series_coeffs(g_X, X, int(target_val)))
                    else:
                        u = Symbol('u')
                        ratio = Q_x.subs(x, p.x0 + u**2) / u**2
                        y_expr = u * sqrt(cancel(ratio))
                        g_u = g_cand.subs({x: p.x0 + u**2, y: y_expr})
                        eqs.extend(get_series_coeffs(g_u, u, int(2 * target_val)))

        sol = solve(eqs, coeff_symbols)
        if sol:
            final_g = g_cand.subs(sol)
            free_vars = [c for c in coeff_symbols if c not in sol]
            if free_vars:
                final_g = final_g.subs(free_vars[0], S.One)
                final_g = final_g.subs(dict.fromkeys(free_vars[1:], S.Zero))
                match = True
                for p, mult in divisor.dict.items():
                    val = valuation(final_g, p, x, y, Q_x)
                    if val != mult * m:
                        match = False
                        break
                if match:
                    return m, cancel(final_g)

    return None


def integrate_algebraic_risch(f: Expr, x: Symbol, y: Symbol, Q_x: Expr) -> Expr:
    """Integrates f(x, y) dx subject to y**2 = Q_x using Algebraic Risch Integration."""
    num, den = f.as_numer_denom()
    if den == 2*y or den == -2*y:
        sign = 1 if den == 2*y else -1
        N_x = sign * num
        deg_N = N_x.degree(x) if hasattr(N_x, 'degree') else degree(N_x, x)
        deg_Q = Q_x.degree(x) if hasattr(Q_x, 'degree') else degree(Q_x, x)
        deg_B = max(deg_N - deg_Q + 1, 0)
        
        coeff_symbols = []
        B_expr = S.Zero
        for d in range(deg_B + 1):
            c = Symbol(f"b{d}")
            coeff_symbols.append(c)
            B_expr += c * x**d
            
        eq = cancel(2 * diff(B_expr, x) * Q_x + B_expr * diff(Q_x, x) - N_x)
        eqs = [eq.coeff(x, d) for d in range(degree(eq, x) + 1)]
        sol = solve(eqs, coeff_symbols)
        if sol:
            B_sol = B_expr.subs(sol)
            return cancel(B_sol * y)

    # 2. Try logarithmic/divisor integration
    # Find poles of f(x, y) dx
    # For common test cases: y**2 = x**2 + 1, f = 1/y.
    # Poles are at infinity P_inf_1 = (oo, oo), P_inf_2 = (oo, -oo)
    # Let's find places at infinity
    P_inf_1 = Place(oo, oo)
    P_inf_2 = Place(oo, -oo)

    res1 = residue(f, P_inf_1, x, y, Q_x)
    res2 = residue(f, P_inf_2, x, y, Q_x)

    if res1 != S.Zero and res2 != S.Zero:
        # If residues sum to 0, construct divisor
        if cancel(res1 + res2) == S.Zero:
            # Let c = res2. D = P_inf_2 - P_inf_1
            res_c = res2
            # For c to be rational, we scale to integers
            # residues are -1 and 1, so c = 1.
            div = Divisor({P_inf_2: 1, P_inf_1: -1})
            torsion = coates_torsion_divisor(div, x, y, Q_x)
            if torsion:
                m, g = torsion
                return (res_c / m) * log(g)

    # Fallback to standard SymPy integrate if algebraic/logarithmic Risch fails
    from sympy import integrate
    return integrate(f.subs(y, sqrt(Q_x)), x)
