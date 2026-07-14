from __future__ import annotations

"""Cylindrical Algebraic Decomposition (CAD) solver.

This module provides functions to compute the Cylindrical Algebraic
Decomposition of a set of polynomials, and to solve arbitrary systems
of polynomial inequalities/equalities over the real numbers.
"""

import math
from sympy import (
    symbols, Poly, discriminant, resultant, S,
    And, Or, Not, Lt, Le, Gt, Ge, Eq, Ne, primitive
)
from sympy.polys.polytools import real_roots

class CADCell:
    """Represents a cell in the Cylindrical Algebraic Decomposition."""
    def __init__(self, sample_point, cell_type):
        self.sample_point = sample_point  # tuple of Rational/Algebraic numbers
        self.cell_type = cell_type        # tuple of 'point' or 'interval'

    def __repr__(self):
        return f"CADCell(sample_point={self.sample_point}, cell_type={self.cell_type})"


def projection(F, var):
    """Computes Collins' projection operator on a set of polynomials F with respect to var."""
    proj_set = set()
    polys_with_var = []
    for f in F:
        if not f.has(var):
            proj_set.add(f)
        else:
            polys_with_var.append(f)

    # 1. Add coefficients and discriminant
    for f in polys_with_var:
        p = Poly(f, var)
        for coeff in p.all_coeffs():
            if coeff != 0 and not coeff.is_number:
                proj_set.add(coeff)
        d = discriminant(f, var)
        if d != 0 and not d.is_number:
            proj_set.add(d)

    # 2. Add pairwise resultants
    n = len(polys_with_var)
    for i in range(n):
        for j in range(i + 1, n):
            r = resultant(polys_with_var[i], polys_with_var[j], var)
            if r != 0 and not r.is_number:
                proj_set.add(r)

    # 3. Clean up and normalize polynomials
    final_set = set()
    for expr in proj_set:
        try:
            _, prim = primitive(expr)
            p_prim = Poly(prim)
            if p_prim.LC() < 0:
                prim = -prim
            final_set.add(prim)
        except Exception:  # noqa: BLE001
            final_set.add(expr)

    return list(final_set)


def get_rational_between(a, b):
    """Finds a simple rational number strictly between two real algebraic numbers a and b."""
    from sympy import nsimplify
    a_f = float(S(a).evalf())
    b_f = float(S(b).evalf())
    mid = (a_f + b_f) / 2
    for tol in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]:
        q = nsimplify(mid, rational=True, tolerance=tol)
        if a_f < float(q) < b_f:
            return q
    return nsimplify(mid, rational=True)


def get_lower_rational(beta):
    """Returns an integer strictly below beta."""
    val = float(S(beta).evalf())
    return S(math.floor(val) - 1)


def get_upper_rational(beta):
    """Returns an integer strictly above beta."""
    val = float(S(beta).evalf())
    return S(math.ceil(val) + 1)


def isolate_algebraic_roots(f, sample_point, variables, var):
    """Isolates real roots of a polynomial f evaluated at a lower-dimensional sample point."""
    subs_dict = {}
    alg_coords = []
    alg_vars = []
    for v, val in zip(variables[:-1], sample_point):
        if val.is_number and val.is_Rational:
            subs_dict[v] = val
        else:
            alg_coords.append(val)
            alg_vars.append(v)

    f_subs = f.subs(subs_dict)

    if not alg_coords:
        return sorted(set(real_roots(f_subs)))

    from sympy.polys.numberfields import primitive_element
    t = symbols('t')
    try:
        minpoly, lincomb, reps = primitive_element(alg_coords, x=t, ex=True)
    except Exception:  # noqa: BLE001
        raise NotImplementedError(f"Could not find primitive element for coordinates: {alg_coords}")

    # Compute minimal polynomial of the primitive element
    M_t = minpoly

    # Map each alg_var to its representation in terms of t
    t_subs = {}
    for v, rep in zip(alg_vars, reps):
        deg = len(rep) - 1
        expr_t = sum(S(c) * t**(deg - j) for j, c in enumerate(rep))
        t_subs[v] = expr_t

    # Substitute representations into f_subs
    G_t_var = f_subs.subs(t_subs)

    # Compute resultant of G_t_var and M_t with respect to t
    R_var = resultant(G_t_var, M_t, t)

    # Find candidate real roots
    cand_roots = real_roots(R_var)

    # Float-evaluate the specific primitive element
    theta_val = sum(c * val for c, val in zip(lincomb, alg_coords))
    theta_num = theta_val.evalf(50)

    # Filter candidate roots using G_t_var(theta, root) == 0
    filtered_roots = []
    for root in cand_roots:
        root_num = root.evalf(50)
        try:
            val = G_t_var.evalf(50, subs={t: theta_num, var: root_num})
            if abs(val) < 1e-40:
                filtered_roots.append(root)
        except Exception:  # noqa: BLE001
            pass

    return sorted(set(filtered_roots))


def cad(system, variables):
    """Computes the Cylindrical Algebraic Decomposition cells for a system of polynomials."""
    # 1. Projection Phase
    F = [system]
    for i in range(len(variables) - 1, 0, -1):
        F.append(projection(F[-1], variables[i]))
    F.reverse()

    # 2. Base Phase (1D)
    F1 = F[0]
    roots_1d = []
    for f in F1:
        if f.has(variables[0]):
            roots_1d.extend(real_roots(f))
    sorted_roots = sorted(set(roots_1d))

    cells = []
    if not sorted_roots:
        cells.append(CADCell(sample_point=(S(0),), cell_type=('interval',)))
    else:
        # First interval
        cells.append(CADCell(sample_point=(get_lower_rational(sorted_roots[0]),), cell_type=('interval',)))
        for i, root in enumerate(sorted_roots):
            cells.append(CADCell(sample_point=(root,), cell_type=('point',)))
            if i < len(sorted_roots) - 1:
                cells.append(CADCell(
                    sample_point=(get_rational_between(root, sorted_roots[i+1]),),
                    cell_type=('interval',)
                ))
        # Last interval
        cells.append(CADCell(sample_point=(get_upper_rational(sorted_roots[-1]),), cell_type=('interval',)))

    # 3. Lifting Phase
    for k in range(1, len(variables)):
        var = variables[k]
        new_cells = []
        for cell in cells:
            s = cell.sample_point
            level_roots = []
            for f in F[k]:
                try:
                    level_roots.extend(isolate_algebraic_roots(f, s, variables[:k+1], var))
                except Exception:  # noqa: BLE001
                    pass
            sorted_level_roots = sorted(set(level_roots))

            if not sorted_level_roots:
                new_cells.append(CADCell(sample_point=s + (S(0),), cell_type=cell.cell_type + ('interval',)))
            else:
                # First interval
                new_cells.append(CADCell(
                    sample_point=s + (get_lower_rational(sorted_level_roots[0]),),
                    cell_type=cell.cell_type + ('interval',)
                ))
                for i, root in enumerate(sorted_level_roots):
                    new_cells.append(CADCell(sample_point=s + (root,), cell_type=cell.cell_type + ('point',)))
                    if i < len(sorted_level_roots) - 1:
                        new_cells.append(CADCell(
                            sample_point=s + (get_rational_between(root, sorted_level_roots[i+1]),),
                            cell_type=cell.cell_type + ('interval',)
                        ))
                # Last interval
                new_cells.append(CADCell(
                    sample_point=s + (get_upper_rational(sorted_level_roots[-1]),),
                    cell_type=cell.cell_type + ('interval',)
                ))
        cells = new_cells

    return cells


def evaluate_relation(rel, subs_dict):
    """Evaluates the truth value of a relational expression under subs_dict."""
    expr = (rel.lhs - rel.rhs).subs(subs_dict)
    if expr.is_number:
        if expr.is_zero:
            sign = 0
        elif expr.is_positive:
            sign = 1
        elif expr.is_negative:
            sign = -1
        else:
            val = expr.evalf(50)
            if abs(val) < 1e-40:
                sign = 0
            elif val > 0:
                sign = 1
            else:
                sign = -1
    else:
        val = expr.evalf(50)
        if abs(val) < 1e-40:
            sign = 0
        elif val > 0:
            sign = 1
        else:
            sign = -1

    if isinstance(rel, Lt):
        return sign < 0
    elif isinstance(rel, Le):
        return sign <= 0
    elif isinstance(rel, Gt):
        return sign > 0
    elif isinstance(rel, Ge):
        return sign >= 0
    elif isinstance(rel, Eq):
        return sign == 0
    elif isinstance(rel, Ne):
        return sign != 0
    return False


def evaluate_formula(formula, subs_dict):
    """Recursively evaluates the truth value of a boolean formula under subs_dict."""
    if isinstance(formula, And):
        return all(evaluate_formula(arg, subs_dict) for arg in formula.args)
    elif isinstance(formula, Or):
        return any(evaluate_formula(arg, subs_dict) for arg in formula.args)
    elif isinstance(formula, Not):
        return not evaluate_formula(formula.args[0], subs_dict)
    return evaluate_relation(formula, subs_dict)


def extract_polys(formula):
    """Extracts all polynomial expressions from a boolean formula of relations."""
    if isinstance(formula, (And, Or, Not)):
        polys = []
        for arg in formula.args:
            polys.extend(extract_polys(arg))
        return list(set(polys))
    return [formula.lhs - formula.rhs]


def solve_cad(formula, variables):
    """Solves a system of inequalities/equalities using CAD.

    Returns the list of CAD cells where the formula is True.
    """
    polys = extract_polys(formula)
    cells = cad(polys, variables)
    matching_cells = []
    for cell in cells:
        subs_dict = dict(zip(variables, cell.sample_point))
        if evaluate_formula(formula, subs_dict):
            matching_cells.append(cell)
    return matching_cells


class Exists:
    """Represents an existential quantifier: Exists(x, formula)."""
    def __init__(self, variable, formula):
        self.variable = variable
        self.formula = formula

    def __repr__(self):
        return f"Exists({self.variable}, {self.formula})"


class ForAll:
    """Represents a universal quantifier: ForAll(x, formula)."""
    def __init__(self, variable, formula):
        self.variable = variable
        self.formula = formula

    def __repr__(self):
        return f"ForAll({self.variable}, {self.formula})"


def parse_quantifiers(formula):
    """Parses nested quantifiers in a formula.

    Returns (quantifiers_list, quantifier_free_formula, free_vars)
    """
    quantifiers_list = []
    curr = formula
    while isinstance(curr, (Exists, ForAll)):
        q_type = 'exists' if isinstance(curr, Exists) else 'forall'
        quantifiers_list.append((q_type, curr.variable))
        curr = curr.formula

    all_vars = sorted(curr.free_symbols, key=lambda s: s.name)
    quantified_vars = [q[1] for q in quantifiers_list]
    free_vars = [v for v in all_vars if v not in quantified_vars]

    return quantifiers_list, curr, free_vars


def simplify_relations(formula):
    """Simplifies logical combinations of relations (e.g. Lt(x, 0) | Eq(x, 0) -> Le(x, 0))."""
    if isinstance(formula, Or):
        args = list(formula.args)
        lts = [a for a in args if isinstance(a, Lt)]
        gts = [a for a in args if isinstance(a, Gt)]
        eqs = [a for a in args if isinstance(a, Eq)]

        new_args = []
        combined = set()

        for eq in eqs:
            matched = False
            for lt in lts:
                if lt.lhs == eq.lhs and lt.rhs == eq.rhs:
                    new_args.append(Le(eq.lhs, eq.rhs))
                    combined.add(eq)
                    combined.add(lt)
                    matched = True
                    break
                elif lt.lhs == eq.rhs and lt.rhs == eq.lhs:
                    new_args.append(Ge(eq.lhs, eq.rhs))
                    combined.add(eq)
                    combined.add(lt)
                    matched = True
                    break
            if not matched:
                for gt in gts:
                    if gt.lhs == eq.lhs and gt.rhs == eq.rhs:
                        new_args.append(Ge(eq.lhs, eq.rhs))
                        combined.add(eq)
                        combined.add(gt)
                        matched = True
                        break
                    elif gt.lhs == eq.rhs and gt.rhs == eq.lhs:
                        new_args.append(Le(eq.lhs, eq.rhs))
                        combined.add(eq)
                        combined.add(gt)
                        matched = True
                        break

        for a in args:
            if a not in combined:
                new_args.append(a)

        if len(new_args) == 1:
            return new_args[0]
        return Or(*new_args)
    return formula


def solve_qe(formula):
    """Performs Quantifier Elimination (QE) using Cylindrical Algebraic Decomposition."""
    quantifiers_list, qf_formula, free_vars = parse_quantifiers(formula)

    if not quantifiers_list:
        return qf_formula

    quantified_vars = [q[1] for q in quantifiers_list]
    variables = free_vars + quantified_vars

    polys = extract_polys(qf_formula)
    cells = cad(polys, variables)

    cell_truth = {}
    for cell in cells:
        subs_dict = dict(zip(variables, cell.sample_point))
        cell_truth[cell.sample_point] = evaluate_formula(qf_formula, subs_dict)

    r = len(variables)
    k = len(free_vars)
    for j in range(r, k, -1):
        q_type, q_var = quantifiers_list[j - k - 1]
        groups = {}
        for s_point, val in cell_truth.items():
            prefix = s_point[:j - 1]
            groups.setdefault(prefix, []).append(val)

        new_truth = {}
        for prefix, vals in groups.items():
            if q_type == 'exists':
                new_truth[prefix] = any(vals)
            else:
                new_truth[prefix] = all(vals)
        cell_truth = new_truth

    if k == 0:
        return S.true if cell_truth.get((), False) else S.false

    F = [polys]
    for i in range(len(variables) - 1, 0, -1):
        F.append(projection(F[-1], variables[i]))
    F.reverse()

    F_k = F[k - 1]

    cells_k = cad(F_k, variables[:k])

    true_descriptions = []
    for cell_k in cells_k:
        s_k = cell_k.sample_point
        if cell_truth.get(s_k, False):
            desc_args = []
            subs_dict = dict(zip(variables[:k], s_k))
            for f in F_k:
                val = f.subs(subs_dict)
                if val.is_number:
                    if val.is_zero:
                        sign = 0
                    elif val.is_positive:
                        sign = 1
                    else:
                        sign = -1
                else:
                    val_num = val.evalf(50)
                    if abs(val_num) < 1e-40:
                        sign = 0
                    elif val_num > 0:
                        sign = 1
                    else:
                        sign = -1

                if sign == 0:
                    desc_args.append(Eq(f, 0))
                elif sign > 0:
                    desc_args.append(Gt(f, 0))
                else:
                    desc_args.append(Lt(f, 0))

            if len(desc_args) == 1:
                true_descriptions.append(desc_args[0])
            elif len(desc_args) > 1:
                true_descriptions.append(And(*desc_args))
            else:
                true_descriptions.append(S.true)

    if not true_descriptions:
        return S.false

    if len(true_descriptions) == 1:
        result_formula = true_descriptions[0]
    else:
        result_formula = Or(*true_descriptions)

    return simplify_relations(result_formula)
