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

def inverse_function(expr, var):
    return sympy.solveset(expr - result, var)

def induction(expr, var):
    return sympy.rsolve(expr.subs(var, expression_iterate(iteration_counter - 1)) - expression_iterate(iteration_counter), expression_iterate(iteration_counter), [var])

def decontextualize_conditions(piecewise):
    pairs = list()
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
        return tuple(result)[0]
    elif isinstance(result, sympy.ImageSet):
        return result.lamda(k_symbol)
    elif isinstance(result, sympy.ConditionSet):
        result = result.condition

        if isinstance(result, sympy.And):
            exprs = result.args
        else:
            exprs = (result,)

        value = 0
        for expr in exprs:
            if isinstance(expr, (sympy.LessThan, sympy.StrictLessThan)):
                pass
            elif isinstance(expr, sympy.GreaterThan):
                value = sympy.ceiling(expr.rhs)
            elif isinstance(expr, sympy.StrictGreaterThan):
                value = sympy.Piecewise((expr.rhs + 1, sympy.Equality(sympy.frac(expr.rhs), 0)), (sympy.ceiling(expr.rhs), sympy.true))
            elif isinstance(expr, sympy.Equality):
                value = expr.rhs
            else:
                raise NotImplementedError

        return value
    elif isinstance(result, (sympy.Interval, sympy.Intersection)):
        if isinstance(result, sympy.Intersection):
            assert tuple(result.args)[0] is sympy.Integers
            result = tuple(result.args)[1]

        return result.left
    else:
        raise NotImplementedError

def _solve_k(expr, var, conditions):
    if isinstance(conditions, sympy.And):
        value = sympy.sympify(0)
        for subcondition in map(partial(_solve_k, expr, var), conditions.args):
            if subcondition.has(k_symbol):
                value, subcondition = subcondition, value

            if value.has(k_symbol):
                inverse = inverse_function(value, k_symbol)
                if not (isinstance(inverse, sympy.FiniteSet) and len(inverse) == 1):
                    raise NotImplementedError
                else:
                    inverse = tuple(inverse)[0]

                inverse = inverse.subs(result, subcondition)
                value = value.subs(k_symbol, sympy.Max(sympy.ceiling(inverse), 0))
            elif value == 0:
                value = subcondition
            else:
                value = sympy.Max(value, subcondition)

        return value.subs(k_symbol, 0)
    elif isinstance(conditions, sympy.Or):
        return sympy.Min(*map(partial(_solve_k, expr, var), conditions.args)).subs(k_symbol, 0)
    else:
        return __solve_k(expr, var, conditions)

def solve_k(expr, var, conditions):
    return sympy.ceiling(_solve_k(expr, var, conditions).subs(k_symbol, 0))

class RPFMode(Enum):
    SIMPLIFY = 1
    SOLVE = 2

def rpf(piecewise, func, mode):
    if len(func.args) > 1:
        return piecewise

    piecewise = sympy.simplify(piecewise)
    func = func.func

    new_pairs = list()
    exit_points = get_exit_points(piecewise, func)

    for pair in decontextualize_conditions(piecewise):
        if pair not in exit_points:
            try:
                var = tuple(pair[0].free_symbols)[0]
                for exit_point in exit_points:
                    expr, after = transform_nonexit_term(pair[0], func)
                    if sympy.ask(sympy.simplify(sympy.Or(pair[1].subs(var, expr), exit_point[1].subs(var, expr))), sympy.simplify(pair[1])) is True:
                        induced_expr = induction(expr, var)
                        k = solve_k(induction(expr, var), var, exit_point[1])
                        new_expr = induced_expr.subs(iteration_counter, k)

                        induced_after = induction(after, result2)
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

    return sympy.Piecewise(*new_pairs, *exit_points)

def rpfsimplify(piecewise, func):
    return rpf(piecewise, func, RPFMode.SIMPLIFY)

def rpfsolve(piecewise, func):
    r"""Converts a piecewise recursive function into a non-recursive expression

    Parameters
    ==========

    piecewise : Expr or a relational.
        The target piecewise
    func : Function
        The name and parameters of the function

    Limitations
    ==========

    It cannot solve expressions with more than one parameter. Piecewise conditions that are recursive must only result in the same piecewise condition being invoked and one non-recursive condition.
    """

    return rpf(piecewise, func, RPFMode.SOLVE)