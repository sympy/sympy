"""Recursive Piecewise Function Solver"""

from functools import partial
from enum import Enum

import sympy

def DummyFunction(name):
    return sympy.Function(sympy.Dummy(name))

result_symbol = sympy.Dummy("result")
result2_symbol = sympy.Dummy("result2_symbol")
expression_iterate = DummyFunction("expression_iterate")
iteration_counter = sympy.Dummy("iteration_counter")
meta_k = sympy.Dummy("meta_k")

x_ = sympy.Dummy("x")
open_ceiling = sympy.Lambda(x_, sympy.floor(x_) + 1)
open_floor = sympy.Lambda(x_, sympy.ceiling(x_) - 1)
del x_

def simplify_piecewise(piecewise):
    expr = sympy.simplify(piecewise)

    if isinstance(expr, sympy.Piecewise):
        return sympy.Piecewise(*((sympy.simplify(term), sympy.simplify(condition)) for term, condition in sympy.simplify(piecewise).args))

    return expr

def inverse_function(expr, var):
    result = sympy.solveset(expr - result_symbol, var)

    if not (result.is_finite and len(result) == 1):
        raise NotImplementedError

    return tuple(result)[0]

def induction(expr, var):
    return sympy.rsolve(expr.subs(var, expression_iterate(iteration_counter - 1)) - expression_iterate(iteration_counter), expression_iterate(iteration_counter), [var])

def univariate_induction(exprs, parameters, var, solved, excluded=None):
    if var in solved:
        return solved[var]

    expr = exprs[parameters.index(var)]

    if excluded is not None:
        for excluded_var in excluded:
            if expr.has(excluded_var):
                raise NotImplementedError("cyclic parameter dependencies are not supported")

    if excluded is None:
        excluded = (var,)

    for var2 in (expr.free_symbols & set(parameters)) - {var}:
        expr = expr.subs(var2, univariate_induction(exprs, parameters, var2, solved, excluded=excluded).subs(iteration_counter, iteration_counter - 1))

    expr = induction(expr, var)

    solved[var] = expr
    return expr

def after_induction(expr, induced_expr_map, k):
    expr = safe_subs(expr.subs(result2_symbol, expression_iterate(iteration_counter - 1)), {symbol: induced_expr.subs(iteration_counter, k - iteration_counter) for symbol, induced_expr in induced_expr_map.items()})

    try:
        return sympy.rsolve(expr.subs(result2_symbol, expression_iterate(iteration_counter - 1)) - expression_iterate(iteration_counter), expression_iterate(iteration_counter), [result2_symbol])
    except (NotImplementedError, ValueError) as error:
        a = sympy.Wild("a")
        b = sympy.Wild("b")

        match = expr.match(expression_iterate(a) + b)
        if match:
            return result2_symbol + sympy.Sum(match[b], (iteration_counter, 1, k))

        match = expr.match(expression_iterate(a) * b)
        if match:
            return result2_symbol + sympy.Product(match[b], (iteration_counter, 1, k))

        raise NotImplementedError from error

def decontextualize_conditions(piecewise):
    pairs = list()
    precondition = sympy.true
    for pair in piecewise.args:
        pairs.append((pair[0], sympy.simplify(sympy.And(pair[1], precondition))))
        precondition = sympy.And(precondition, sympy.Not(pair[1]))

    return tuple(pairs)

def get_exit_points(piecewise_args, func):
    return tuple(filter(lambda pair: not pair[0].has(func), piecewise_args))

def transform_nonexit_term(expr, func):
    function_exprs = tuple(expr.atoms(func))
    function_expr = function_exprs[0]
    if not all(map(function_expr.equals, function_exprs[1:])):
        raise NotImplementedError("recursive terms cannot branch to more than one recursive term or have the possibility of doing so")

    return function_expr.args, expr.replace(function_expr, result2_symbol)

def safe_subs(expr, replacements):
    symbols_replaced = tuple(expr.free_symbols & set(replacements.keys()))

    temp_symbols = list()
    for symbol in symbols_replaced:
        temp_symbol = sympy.Dummy()
        expr = expr.subs(symbol, temp_symbol)
        temp_symbols.append(temp_symbol)

    for index, temp_symbol in enumerate(temp_symbols):
        expr = expr.subs(temp_symbol, replacements[symbols_replaced[index]])

    return expr

def __solve_k_process_remove_integer_intersections(result):
    if isinstance(result, sympy.Intersection):
        result = list(result.args)
        if sympy.S.Integers in result:
            result.remove(sympy.S.Integers)
            if len(result) == 1:
                return result[0]
            else:
                raise NotImplementedError
        else:
            raise NotImplementedError

    return result

class KConstraintType(Enum):
    EQUAL = 1
    GREATER = 2
    LESSER = 3
    GENERATOR = 4
    PRECONDITION = 5

def __solve_k_process(result):
    if isinstance(result, sympy.FiniteSet):
        return [[(item, KConstraintType.EQUAL)] for item in result]
    elif isinstance(result, sympy.ImageSet):
        return [[(result.lamda(meta_k), KConstraintType.GENERATOR)]]
    elif isinstance(result, sympy.Interval):
        return [[(open_ceiling(result.left) if result.left_open else sympy.ceiling(result.left), KConstraintType.GREATER), (open_floor(result.right) if result.right_open else sympy.floor(result.right), KConstraintType.LESSER)]]
    elif isinstance(result, sympy.Intersection):
        return __solve_k_process(__solve_k_process_remove_integer_intersections(result))
    elif isinstance(result, sympy.Union):
        return sum(map(__solve_k_process, result.args), start=list())
    else:
        raise NotImplementedError

def __solve_k(induced_expr_map, condition):
    condition = safe_subs(condition, induced_expr_map)

    try:
        return __solve_k_process(sympy.solveset(condition, iteration_counter, domain=sympy.S.Integers))
    except NotImplementedError as error:
        if isinstance(condition, (sympy.core.relational.GreaterThan, sympy.core.relational.LessThan, sympy.core.relational.StrictGreaterThan, sympy.core.relational.StrictLessThan)):
            root_expr = condition.lhs - condition.rhs
            roots = __solve_k_process_remove_integer_intersections(sympy.solveset(root_expr, iteration_counter, domain=sympy.S.Integers))

            if isinstance(roots, sympy.FiniteSet):
                out = list()
                for root in roots + {-sympy.oo}:
                    roots2 = (roots - {root}) + {sympy.oo}
                    min_root2 = sympy.Min(*(sympy.Piecewise((root2, sympy.GreaterThan(root2, root)), (sympy.oo, sympy.true)) for root2 in roots2))

                    for root2 in roots2:
                        if {root, root2} == {-sympy.oo, sympy.oo}:
                            continue

                        out.append([(open_ceiling(root) if isinstance(condition, (sympy.core.relational.StrictGreaterThan, sympy.core.relational.StrictLessThan)) else sympy.ceiling(root), KConstraintType.GREATER), (open_floor(root2) if isinstance(condition, (sympy.core.relational.StrictGreaterThan, sympy.core.relational.StrictLessThan)) else sympy.floor(root2), KConstraintType.LESSER), (sympy.And(sympy.Equality(root2, min_root2), type(condition)(root_expr.subs(iteration_counter, (root + root2) / 2), 0)), KConstraintType.PRECONDITION)])

                return out

        raise NotImplementedError from error

def _solve_k(induced_expr_map, conditions):
    if isinstance(conditions, sympy.And):
        results = tuple(map(partial(_solve_k, induced_expr_map), conditions.args))

        permutations = [list()]
        for item in results:
            new_permutations = list()
            for item2 in permutations:
                for item3 in item:
                    new_permutations.append(item2 + item3)
            permutations = new_permutations

        return permutations
    elif isinstance(conditions, sympy.Or):
        return sum(map(partial(_solve_k, induced_expr_map), conditions.args), start=list())
    else:
        return __solve_k(induced_expr_map, conditions)

def solve_k(induced_expr_map, conditions):
    permutations = _solve_k(induced_expr_map, conditions)

    out = list()
    for permutation in permutations:
        values, types = zip(*permutation)

        precondition = sympy.true

        if KConstraintType.EQUAL in types:
            expr = values[types.index(KConstraintType.EQUAL)]

            for index2, (value, type_) in enumerate(permutation):
                if type_ is KConstraintType.EQUAL:
                    precondition = sympy.And(precondition, sympy.Equality(sympy.frac(value), 0))
                    precondition = sympy.And(precondition, sympy.Equality(value, expr))
                    precondition = sympy.And(precondition, sympy.GreaterThan(value, 1))
                elif type_ is KConstraintType.GREATER:
                    precondition = sympy.And(precondition, sympy.GreaterThan(expr, value))
                elif type_ is KConstraintType.LESSER:
                    precondition = sympy.And(precondition, sympy.LessThan(expr, value))
                elif type_ is KConstraintType.GENERATOR:
                    precondition = sympy.And(precondition, sympy.Equality(sympy.frac(inverse_function(value, meta_k).subs(result_symbol, expr)), 0))
                    precondition = sympy.And(precondition, sympy.GreaterThan(inverse_function(value, meta_k).subs(result_symbol, expr), 0))
                elif type_ is KConstraintType.PRECONDITION:
                    precondition = sympy.And(precondition, value)
        else:
            if types.count(KConstraintType.GENERATOR) > 1:
                raise NotImplementedError

            greater_values, _ = zip(*filter(lambda pair: pair[1] is KConstraintType.GREATER, permutation))
            if greater_values:
                expr = sympy.Max(1, *greater_values)
            else:
                expr = sympy.sympify(1)

            if KConstraintType.GENERATOR in types:
                generator = values[types.index(KConstraintType.GENERATOR)]
                expr = generator.subs(meta_k, sympy.Max(sympy.ceiling(inverse_function(generator, meta_k).subs(result_symbol, expr)), 0))

            for value, _ in filter(lambda pair: pair[1] is KConstraintType.LESSER, permutation):
                precondition = sympy.And(precondition, sympy.LessThan(expr, value))

            for value, _ in filter(lambda pair: pair[1] is KConstraintType.PRECONDITION, permutation):
                precondition = sympy.And(precondition, value)

        out.append((expr, precondition))

    return out

class RPFMode(Enum):
    SIMPLIFY = 1
    SOLVE = 2

def rpf(piecewise, func, mode):
    parameters = func.args
    func = func.func

    piecewise_args = decontextualize_conditions(piecewise)

    simplifed = True
    while simplifed:
        simplifed = False

        exit_points = get_exit_points(piecewise_args, func)

        new_pairs = list()
        for pair in piecewise_args:
            if pair not in exit_points:
                try:
                    exprs, after = transform_nonexit_term(pair[0], func)
                    next_expr_map = {var: exprs[index] for index, var in enumerate(parameters)}

                    solved = dict()
                    for var in parameters:
                        univariate_induction(exprs, parameters, var, solved)

                    possible_exit_points = list()
                    exit_point_conditions = list()
                    for exit_point in exit_points:
                        if sympy.ask(sympy.simplify(safe_subs(exit_point[1], next_expr_map))) is not False:
                            possible_exit_points.append(exit_point)
                            exit_point_conditions.append(safe_subs(exit_point[1], next_expr_map))

                    if sympy.ask(sympy.simplify(sympy.Or(safe_subs(pair[1], next_expr_map), *exit_point_conditions)), sympy.simplify(pair[1])) is True:
                        new_terms_and_ks = list()
                        checked_ks = list()

                        for exit_point in possible_exit_points:
                            for k, k_precondition in solve_k(solved, exit_point[1]):
                                condition2 = sympy.And(safe_subs(exit_point[1], solved).subs(iteration_counter, k), k_precondition)

                                induced_after = after_induction(after, solved, k)
                                new_terms_and_ks.append(((induced_after.subs(result2_symbol, safe_subs(exit_point[0], solved).subs(iteration_counter, k)).subs(iteration_counter, k), sympy.And(pair[1], condition2)), k))
                                checked_ks.append(sympy.Piecewise((k, condition2), (sympy.oo, sympy.true)))

                        min_k = sympy.Min(*checked_ks)

                        for (term, condition), k in new_terms_and_ks:
                            new_pairs.append((term, sympy.And(condition, sympy.Equality(k, min_k))))

                        simplifed = True
                    else:
                        if mode is RPFMode.SIMPLIFY:
                            new_pairs.append(pair)
                        else:
                            raise NotImplementedError("recursive terms cannot branch to more than one recursive term or have the possibility of doing so")
                except NotImplementedError as error:
                    if mode is RPFMode.SIMPLIFY:
                        new_pairs.append(pair)
                    else:
                        raise NotImplementedError from error
            else:
                new_pairs.append(pair)

        piecewise_args = new_pairs

    return new_pairs

def rpfsimplify(piecewise, func):
    r"""Attempts to simplify recursive terms in piecewise functions into a non-recursive ones

    Parameters
    ==========

    piecewise : Piecewise
        The target piecewise
    func : Function
        The name and parameters of the function

    Returns
    =======

    Piecewise
        The simplified piecewise.

    Limitations
    ======

    It cannot solve expressions with the cyclically dependent parameters (e.g. f(x + 1, y + 1) and f(x + y, y + 1) are allowed, but f(x + y, y + x) is not). Piecewise conditions that are recursive must not branch (i.e. they cannot have the possibility of invoking more than one recursive term [including itself]).

    See Also
    ========

    sympy.solvers.rpfsolve: RPF solver; only returns if all recursive terms can be solved

    """

    return rpf(piecewise, func, RPFMode.SIMPLIFY)

def rpfsolve(piecewise, func):
    r"""Converts a piecewise recursive function into a non-recursive expression

    Parameters
    ==========

    piecewise : Piecewise
        The target piecewise
    func : Function
        The name and parameters of the function

    Returns
    =======

    Piecewise
        A transformed non-recursive version of the input piecewise

    Raises
    ======

    NotImplementedError
        If any of the recursive terms cannot be solved.

    Limitations
    ======

    It cannot solve expressions with the cyclically dependent parameters (e.g. f(x + 1, y + 1) and f(x + y, y + 1) are allowed, but f(x + y, y + x) is not). Piecewise conditions that are recursive must not branch (i.e. they cannot have the possibility of invoking more than one recursive term [including itself]).

    See Also
    ========

    sympy.rpfsimplify: simplifier using rpfsolve; returns even if all terms can't be solved

    """

    return rpf(piecewise, func, RPFMode.SOLVE)
