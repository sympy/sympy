"""Recursive Piecewise Function Solver"""

import sympy

from functools import partial
from enum import Enum

def DummyFunction(name):
    return sympy.Function(sympy.Dummy(name))

result = sympy.Dummy("result")
result2 = sympy.Dummy("result2")
expression_iterate = DummyFunction("expression_iterate")
iteration_counter = sympy.Dummy("iteration_counter")
meta_k = sympy.Dummy("meta_k")

x = sympy.Dummy("x")
open_ceiling_left = sympy.Lambda(x, sympy.floor(x) + 1)
open_ceiling_right = sympy.Lambda(x, sympy.ceiling(x) - 1)

def simplify_piecewise(piecewise):
    expr = sympy.simplify(piecewise)

    if isinstance(expr, sympy.Piecewise):
        return sympy.Piecewise(*((sympy.simplify(term), sympy.simplify(condition)) for term, condition in sympy.simplify(piecewise).args))

    return expr

def inverse_function(expr, var):
    inverse = sympy.solveset(expr - result, var)
    if not (isinstance(inverse, sympy.FiniteSet) and len(inverse) == 1):
        raise NotImplementedError

    return tuple(inverse)[0]

def induction(expr, var):
    return sympy.rsolve(expr.subs(var, expression_iterate(iteration_counter - 1)) - expression_iterate(iteration_counter), expression_iterate(iteration_counter), [var])

def after_induction(expr, induced_expr, k, result_var, func_var):
    expr = expr.subs(result_var, expression_iterate(iteration_counter - 1)).subs(func_var, induced_expr.subs(iteration_counter, k - iteration_counter))

    try:
        return sympy.rsolve(expr.subs(result_var, expression_iterate(iteration_counter - 1)).subs(func_var, induced_expr.subs(iteration_counter, k - iteration_counter)) - expression_iterate(iteration_counter), expression_iterate(iteration_counter), [result_var])
    except (NotImplementedError, ValueError):
        a = sympy.Wild("a")
        b = sympy.Wild("b")

        match = expr.match(expression_iterate(a) + b)
        if match:
            return sympy.Sum(match[b], (iteration_counter, 1, k))

        match = expr.match(expression_iterate(a) * b)
        if match:
            return sympy.Product(match[b], (iteration_counter, 1, k))

        raise NotImplementedError

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
    assert all(map(function_expr.equals, function_exprs[1:]))

    return function_expr.args[0], expr.replace(function_expr, result2)

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

    raise NotImplementedError

class KConstraintType(Enum):
    EQUAL = 1
    GREATER = 2
    LESSER = 3
    GENERATOR = 4

def __solve_k_process(result):
    if isinstance(result, sympy.FiniteSet):
        return [[(item, KConstraintType.EQUAL)] for item in result]
    elif isinstance(result, sympy.ImageSet):
        return [[(result.lamda(meta_k), KConstraintType.GENERATOR)]]
    elif isinstance(result, sympy.Interval):
        return [[(open_ceiling_left(result.left) if result.left_open else sympy.ceiling(result.left), KConstraintType.GREATER), (open_ceiling_right(result.right) if result.right_open else sympy.ceiling(result.right), KConstraintType.LESSER)]]
    elif isinstance(result, sympy.Intersection):
        return __solve_k_process(__solve_k_process_remove_integer_intersections(result))
    elif isinstance(result, sympy.Union):
        return sum(map(__solve_k_process, result.args), start=list())
    else:
        raise NotImplementedError

def __solve_k(expr, var, condition):
    return __solve_k_process(sympy.solveset(condition.subs(var, expr), iteration_counter, domain=sympy.S.Integers))

def _solve_k(expr, var, conditions):
    if isinstance(conditions, sympy.And):
        results = tuple(map(partial(_solve_k, expr, var), conditions.args))

        permutations = [list()]
        for item in results:
            new_permutations = list()
            for item2 in permutations:
                for item3 in item:
                    new_permutations.append(item2 + item3)
            permutations = new_permutations

        return permutations
    elif isinstance(conditions, sympy.Or):
        return sum(map(partial(_solve_k, expr, var), conditions.args), start=list())
    else:
        return __solve_k(expr, var, conditions)

def solve_k(expr, var, conditions):
    permutations = _solve_k(expr, var, conditions)

    out = list()
    for permutation in permutations:
        values, types = zip(*permutation)

        prior_condition = sympy.true

        if KConstraintType.EQUAL in types:
            expr = values[types.index(KConstraintType.EQUAL)]
            for index2, (value, type_) in enumerate(permutation):
                if type_ is KConstraintType.EQUAL:
                    prior_condition = sympy.And(prior_condition, sympy.Equality(sympy.frac(value), 0))
                    prior_condition = sympy.And(prior_condition, sympy.Equality(value, expr))
                elif type_ is KConstraintType.GREATER:
                    prior_condition = sympy.And(prior_condition, sympy.GreaterThan(expr, value))
                elif type_ is KConstraintType.LESSER:
                    prior_condition = sympy.And(prior_condition, sympy.LessThan(expr, value))
                elif type_ is KConstraintType.GENERATOR:
                    prior_condition = sympy.And(prior_condition, sympy.Equality(sympy.frac(inverse_function(value, meta_k).subs(result, expr)), 0))
                    prior_condition = sympy.And(prior_condition, sympy.GreaterThan(inverse_function(value, meta_k).subs(result, expr), 0))
        else:
            if types.count(KConstraintType.GENERATOR) > 1:
                raise NotImplementedError

            greater_values, _ = zip(*filter(lambda pair: pair[1] is KConstraintType.GREATER, permutation))
            if greater_values:
                expr = sympy.Max(*greater_values)
            else:
                expr = sympy.sympify(0)

            for value, _ in filter(lambda pair: pair[1] is KConstraintType.LESSER, permutation):
                prior_condition = sympy.And(prior_condition, sympy.LessThan(expr, value))

            if KConstraintType.GENERATOR in types:
                generator = values[types.index(KConstraintType.GENERATOR)]
                expr = generator.subs(meta_k, sympy.Max(sympy.ceiling(inverse_function(generator, meta_k).subs(result, expr)), 0))

        out.append((expr, prior_condition))

    return out

class RPFMode(Enum):
    SIMPLIFY = 1
    SOLVE = 2

def rpf(piecewise, func, mode):
    if len(func.args) > 1:
        if mode is RPFMode.SIMPLIFY:
            return piecewise
        else:
            raise NotImplementedError

    piecewise = simplify_piecewise(piecewise)
    var = func.args[0]
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
                    expr, after = transform_nonexit_term(pair[0], func)
                    induced_expr = induction(expr, var)
    
                    possible_exit_points = list()
                    exit_point_conditions = list()
                    for exit_point in exit_points:
                        if sympy.ask(sympy.simplify(exit_point[1].subs(var, expr))) is not False:
                            possible_exit_points.append(exit_point)
                            exit_point_conditions.append(exit_point[1].subs(var, expr))
    
                    if sympy.ask(sympy.simplify(sympy.Or(pair[1].subs(var, expr), *exit_point_conditions)), sympy.simplify(pair[1])) is True:
                        new_terms_and_ks = list()
                        checked_ks = list()
    
                        for exit_point in possible_exit_points:
                            for k, k_condition in solve_k(induced_expr, var, exit_point[1]):
                                new_expr = induced_expr.subs(iteration_counter, k)
                                condition2 = sympy.And(exit_point[1].subs(var, new_expr), k_condition, sympy.GreaterThan(k, 0))
    
                                induced_after = after_induction(after, induced_expr, k, result2, var)
                                new_terms_and_ks.append(((induced_after.subs(result2, exit_point[0].subs(var, new_expr)).subs(iteration_counter, k), sympy.And(pair[1], condition2)), k))
                                checked_ks.append(sympy.Piecewise((k, condition2), (sympy.oo, sympy.true)))
    
                        min_k = sympy.Min(*checked_ks)
    
                        for (term, condition), k in new_terms_and_ks:
                            new_pairs.append((term, sympy.And(condition, sympy.Equality(k, min_k))))

                        simplifed = True
                    else:
                        if mode is RPFMode.SIMPLIFY:
                            new_pairs.append(pair)
                        else:
                            raise NotImplementedError
                except NotImplementedError as error:
                    if mode is RPFMode.SIMPLIFY:
                        new_pairs.append(pair)
                    else:
                        raise NotImplementedError from error
            else:
                new_pairs.append(pair)

        piecewise_args = new_pairs

    return simplify_piecewise(sympy.Piecewise(*new_pairs))

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

    It cannot solve expressions with more than one parameter. Piecewise conditions that are recursive must not branch (i.e. they cannot have the possibility of invoking more than one recursive term [including itself]).

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

    It cannot solve expressions with more than one parameter. Piecewise conditions that are recursive must not branch (i.e. they cannot have the possibility of invoking more than one recursive term [including itself]).

    See Also
    ========

    sympy.rpfsimplify: simplifier using rpfsolve; returns even if all terms can't be solved

    """

    return rpf(piecewise, func, RPFMode.SOLVE)
