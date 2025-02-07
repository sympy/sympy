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
k_symbol = sympy.Dummy("k")

x = sympy.Dummy("x")
open_ceiling = sympy.Lambda(x, sympy.Piecewise((x + 1, sympy.Equality(sympy.frac(x), 0)), (sympy.ceiling(x), sympy.true)))

def inverse_function(expr, var):
    return sympy.solveset(expr - result, var)

def induction(expr, var):
    return sympy.rsolve(expr.subs(var, expression_iterate(iteration_counter - 1)) - expression_iterate(iteration_counter), expression_iterate(iteration_counter), [var])

def after_induction(expr, induced_expr, k, result_var, func_var):
    expr = expr.subs(result_var, expression_iterate(iteration_counter - 1)).subs(func_var, induced_expr.subs(iteration_counter, k - iteration_counter))

    try:
        return sympy.rsolve(expr.subs(result_var, expression_iterate(iteration_counter - 1)).subs(func_var, induced_expr.subs(iteration_counter, k - iteration_counter)) - expression_iterate(iteration_counter), expression_iterate(iteration_counter), [result_var])
    except (NotImplementedError, ValueError):
        x = sympy.Wild("x")
        y = sympy.Wild("y")

        match = expr.match(expression_iterate(x) + y)
        if match:
            return sympy.Sum(match[y], (iteration_counter, 1, k))

        match = expr.match(expression_iterate(x) * y)
        if match:
            return sympy.Product(match[y], (iteration_counter, 1, k))

        raise NotImplementedError

def decontextualize_conditions(piecewise):
    pairs = []
    precondition = sympy.true
    for pair in piecewise.args:
        pairs.append((pair[0], sympy.simplify(sympy.And(pair[1], precondition))))
        precondition = sympy.And(precondition, sympy.Not(pair[1]))

    return tuple(pairs)

def get_exit_points(piecewise, func):
    return tuple(filter(lambda pair: not pair[0].has(func), decontextualize_conditions(piecewise)))

def transform_nonexit_term(expr, func):
    function_exprs = tuple(filter(lambda subexpr: subexpr.func == func, expr.atoms(sympy.Function)))
    function_expr = function_exprs[0]
    assert all(map(function_expr.equals, function_exprs[1:]))

    return function_expr.args[0], expr.replace(function_expr, result2)

def __solve_k(expr, var, condition):
    result = sympy.solveset(condition.subs(var, expr), iteration_counter, domain=sympy.Integers)
    if isinstance(result, sympy.FiniteSet):
        assert len(result) == 1
        return tuple(result)[0], sympy.ceiling
    elif isinstance(result, sympy.ImageSet):
        return result.lamda(k_symbol), sympy.ceiling
    elif isinstance(result, sympy.ConditionSet):
        result = result.condition
        ceiling_func = sympy.ceiling

        if isinstance(result, sympy.And):
            exprs = result.args
        else:
            exprs = (result,)

        value = 0
        for expr in exprs:
            if isinstance(expr, (sympy.LessThan, sympy.StrictLessThan)):
                pass
            elif isinstance(expr, sympy.GreaterThan):
                value = expr.rhs
                ceiling_func = sympy.ceiling
            elif isinstance(expr, sympy.StrictGreaterThan):
                value = expr.rhs
                ceiling_func = open_ceiling
            elif isinstance(expr, sympy.Equality):
                value = expr.rhs
            else:
                raise NotImplementedError

        return value, ceiling_func
    elif isinstance(result, (sympy.Interval, sympy.Intersection)):
        if isinstance(result, sympy.Intersection):
            result = list(result.args)
            if sympy.Integers in result:
                result.remove(sympy.Integers)
                if len(result) == 1:
                    result = result[0]
                else:
                    raise NotImplementedError
            else:
                raise NotImplementedError

        return result.left, open_ceiling if result.left_open else sympy.ceiling
    else:
        raise NotImplementedError

def _solve_k(expr, var, conditions):
    if isinstance(conditions, sympy.And):
        value = sympy.sympify(0)
        previous_ceiling_func = sympy.ceiling
        for index, (subcondition, ceiling_func) in enumerate(map(partial(_solve_k, expr, var), conditions.args)):
            if subcondition.has(k_symbol):
                value, subcondition = subcondition, value
                ceiling_func, previous_ceiling_func = previous_ceiling_func, ceiling_func

            if value.has(k_symbol):
                inverse = inverse_function(value, k_symbol)
                if not (isinstance(inverse, sympy.FiniteSet) and len(inverse) == 1):
                    raise NotImplementedError
                else:
                    inverse = tuple(inverse)[0]

                inverse = inverse.subs(result, subcondition)
                value = value.subs(k_symbol, sympy.Max(ceiling_func(inverse), 0))
            elif value == 0:
                value = subcondition
            else:
                value = sympy.Max(value, subcondition)
                previous_ceiling_func = previous_ceiling_func(ceiling_func)

            if index != len(conditions.args) - 1:
                previous_ceiling_func = ceiling_func

        return value.subs(k_symbol, 0), previous_ceiling_func
    elif isinstance(conditions, sympy.Or):
        return sympy.Min(*(ceiling_func(expr.subs(k_symbol, 0)) for expr, ceiling_func in map(partial(_solve_k, expr, var), conditions.args))), sympy.ceiling
    else:
        return __solve_k(expr, var, conditions)

def solve_k(expr, var, conditions):
    expr, ceiling_func = _solve_k(expr, var, conditions)
    return ceiling_func(expr.subs(k_symbol, 0))

class RPFMode(Enum):
    SIMPLIFY = 1
    SOLVE = 2

def rpf(piecewise, func, mode):
    if len(func.args) > 1:
        if mode is RPFMode.SIMPLIFY:
            return piecewise
        else:
            raise NotImplementedError

    piecewise = sympy.simplify(piecewise)
    var = func.args[0]
    func = func.func

    new_pairs = []
    exit_points = get_exit_points(piecewise, func)

    for pair in decontextualize_conditions(piecewise):
        if pair not in exit_points:
            try:
                for exit_point in exit_points:
                    expr, after = transform_nonexit_term(pair[0], func)
                    if sympy.ask(sympy.simplify(sympy.Or(pair[1].subs(var, expr), exit_point[1].subs(var, expr))), sympy.simplify(pair[1])) is True:
                        induced_expr = induction(expr, var)
                        k = solve_k(induction(expr, var), var, exit_point[1])
                        new_expr = induced_expr.subs(iteration_counter, k)

                        induced_after = after_induction(after, induced_expr, k, result2, var)
                        new_pairs.append((induced_after.subs(result2, exit_point[0].subs(var, new_expr)).subs(iteration_counter, k), sympy.And(pair[1], exit_point[1].subs(var, new_expr))))
                        break
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

    return sympy.Piecewise(*new_pairs)

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

    It cannot solve expressions with more than one parameter. Piecewise conditions that are recursive must only result in the same piecewise condition being invoked and one non-recursive condition.

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

    It cannot solve expressions with more than one parameter. Piecewise conditions that are recursive must only result in the same piecewise condition being invoked and one non-recursive condition.

    See Also
    ========

    sympy.rpfsimplify: simplifier using rpfsolve; returns even if all terms can't be solved

    """

    return rpf(piecewise, func, RPFMode.SOLVE)
