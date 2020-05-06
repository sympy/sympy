#
# The main tests for the code in single.py are currently located in
# sympy/solvers/tests/test_ode.py
#
r"""
This File contains test functions for the individual hints used for solving ODEs.

Examples of each solver will be returned by _get_examples_ode_sol_name_of_solver.

Examples should have a key 'XFAIL' which stores the list of hints if they are
expected to fail for that hint.

Functions that are for internal use:

1) _ode_solver_test(ode_examples) - It takes dictionary of examples returned by
   _get_examples method and tests them with their respective hints.

2) _test_particular_example(our_hint, example_name) - It tests the ODE example corresponding
   to the hint provided.

3) _test_all_hints(runxfail=False) - It is used to test all the examples with all the hints
  currently implemented. It calls _test_all_examples_for_one_hint() which outputs whether the
  given hint functions properly if it classifies the ODE example.
  If runxfail flag is set to True then it will only test the examples which are expected to fail.

  Everytime the ODE of partiular solver are added then _test_all_hints() is to execuetd to find
  the possible failures of different solver hints.

4) _test_all_examples_for_one_hint(our_hint, all_examples) - It takes hint as argument and checks
   this hint against all the ODE examples and gives output as the number of ODEs matched, number
   of ODEs which were solved correctly, list of ODEs which gives incorrect solution and list of
   ODEs which raises exception.

"""
from sympy import (cos, Derivative, Dummy, diff,
    Eq, exp, I, log, pi, Piecewise, Rational, sin, tan,
    sqrt, symbols, Ei, erfi)

from sympy.core import Function, Symbol
from sympy.functions import airyai, airybi, besselj, bessely
from sympy.integrals.risch import NonElementaryIntegral
from sympy.solvers.ode import classify_ode, dsolve
from sympy.solvers.ode.ode import allhints, _remove_redundant_solutions
from sympy.solvers.ode.single import (FirstLinear, ODEMatchError,
    SingleODEProblem, SingleODESolver)

from sympy.solvers.ode.subscheck import checkodesol

from sympy.testing.pytest import raises
import traceback


x = Symbol('x')
y = Symbol('y')
f = Function('f')
g = Function('g')
C1, C2, C3 = symbols('C1:4')


hint_message = """\
Hint did not match the example {example}.

The ODE is:
{eq}.

The expected hint was
{our_hint}\
"""

expected_sol_message = """\
Different solution found from dsolve for example {example}.

The ODE is:
{eq}

The expected solution was
{sol}

What dsolve returned is:
{dsolve_sol}\
"""

checkodesol_msg = """\
solution found is not correct for example {example}.

The ODE is:
{eq}\
"""

dsol_incorrect_msg = """\
solution returned by dsolve is incorrect when using {hint}.

The ODE is:
{eq}

The expected solution was
{sol}

what dsolve returned is:
{dsolve_sol}

You can test this with:

eq = {eq}
sol = dsolve(eq, hint='{hint}')
print(sol)
print(checkodesol(eq, sol))

"""

exception_msg = """\
dsolve raised exception : {e}

when using {hint} for the example {example}

You can test this with:

from sympy.solvers.ode.tests.test_single import _test_an_example

_test_an_example('{hint}', example_name = '{example}')

The ODE is:
{eq}

\
"""

check_hint_msg = """\
Tested hint was : {hint}

Total of {matched} examples matched with this hint.

Out of which {solve} gave correct results.

Examples which gave incorrect results are {unsolve}.

Examples which raised exceptions are {exceptions}
\
"""


def _ode_solver_test(ode_examples):
    our_hint = ode_examples['hint']
    for example in ode_examples['examples']:
        temp = {
            'eq': ode_examples['examples'][example]['eq'],
            'sol': ode_examples['examples'][example]['sol'],
            'XFAIL': ode_examples['examples'][example].get('XFAIL', []),
            'func': ode_examples['examples'][example].get('func',ode_examples['func']),
            'example_name': example,
        }

        result = _test_particular_example(our_hint, temp, solver_flag=True)
        if result['xpass_msg'] != "":
            print(result['xpass_msg'])


def _test_all_hints(runxfail=False):
    all_hints = list(allhints)+["default"]
    all_examples = _get_all_examples()

    for our_hint in all_hints:
        if our_hint.endswith('_Integral') or 'series' in our_hint:
            continue
        _test_all_examples_for_one_hint(our_hint, all_examples, runxfail)


def _test_dummy_sol(expected_sol,dsolve_sol):
    if type(dsolve_sol)==list:
        return any(expected_sol.dummy_eq(sub_dsol) for sub_dsol in dsolve_sol)
    else:
        return expected_sol.dummy_eq(dsolve_sol)


def _test_an_example(our_hint, example_name):
    all_examples = _get_all_examples()
    for example in all_examples:
        if example['example_name'] == example_name:
            _test_particular_example(our_hint, example)


def _test_particular_example(our_hint, ode_example, solver_flag=False):
    eq = ode_example['eq']
    expected_sol = ode_example['sol']
    example = ode_example['example_name']
    xfail = our_hint in ode_example['XFAIL']
    func = ode_example['func']
    result = {'msg': '', 'xpass_msg': ''}
    xpass = True
    if solver_flag:
        if our_hint not in classify_ode(eq, func):
            message = hint_message.format(example=example, eq=eq, our_hint=our_hint)
            raise AssertionError(message)

    if our_hint in classify_ode(eq, func):
        result['match_list'] = example
        try:
            dsolve_sol = dsolve(eq, func, hint=our_hint)

        except Exception as e:
            dsolve_sol = []
            result['exception_list'] = example
            if not solver_flag:
                traceback.print_exc()
                result['msg'] = exception_msg.format(e=str(e), hint=our_hint, example=example, eq=eq)
            xpass = False

        if solver_flag and dsolve_sol!=[]:
            expect_sol_check = False
            if type(dsolve_sol)==list:
                for sub_sol in expected_sol:
                    if sub_sol.has(Dummy):
                        expect_sol_check = not _test_dummy_sol(sub_sol, dsolve_sol)
                    else:
                        expect_sol_check = sub_sol not in dsolve_sol
                    if expect_sol_check:
                        break
            else:
                expect_sol_check = dsolve_sol not in expected_sol
                for sub_sol in expected_sol:
                    if sub_sol.has(Dummy):
                        expect_sol_check = not _test_dummy_sol(sub_sol, dsolve_sol)

            if expect_sol_check:
                message = expected_sol_message.format(example=example, eq=eq, sol=expected_sol, dsolve_sol=dsolve_sol)
                raise AssertionError(message)

            expected_checkodesol = [(True, 0) for i in range(len(expected_sol))]
            if len(expected_sol) == 1:
                expected_checkodesol = (True, 0)

            if checkodesol(eq, dsolve_sol) != expected_checkodesol:
                result['unsolve_list'] = example
                xpass = False
                message = dsol_incorrect_msg.format(hint=our_hint, eq=eq, sol=expected_sol,dsolve_sol=dsolve_sol)
                if solver_flag:
                    message = checkodesol_msg.format(example=example, eq=eq)
                    raise AssertionError(message)
                else:
                    result['msg'] = 'AssertionError: ' + message

        if xpass and xfail:
            result['xpass_msg'] = example + "is now passing for the hint" + our_hint
    return result


def _test_all_examples_for_one_hint(our_hint, all_examples=[], runxfail=None):
    if all_examples == []:
        all_examples = _get_all_examples()
    match_list, unsolve_list, exception_list = [], [], []
    for ode_example in all_examples:
        xfail = our_hint in ode_example['XFAIL']
        if runxfail and not xfail:
            continue
        if xfail:
            continue
        result = _test_particular_example(our_hint, ode_example)
        match_list += result.get('match_list',[])
        unsolve_list += result.get('unsolve_list',[])
        exception_list += result.get('exception_list',[])
        if runxfail is not None:
            msg = result['msg']
            if msg!='':
                print(result['msg'])
            # print(result.get('xpass_msg',''))
    if runxfail is None:
        match_count = len(match_list)
        solved = len(match_list)-len(unsolve_list)-len(exception_list)
        msg = check_hint_msg.format(hint=our_hint, matched=match_count, solve=solved, unsolve=unsolve_list, exceptions=exception_list)
        print(msg)


def test_SingleODESolver():
    # Test that not implemented methods give NotImplementedError
    # Subclasses should override these methods.
    problem = SingleODEProblem(f(x).diff(x), f(x), x)
    solver = SingleODESolver(problem)
    raises(NotImplementedError, lambda: solver.matches())
    raises(NotImplementedError, lambda: solver.get_general_solution())
    raises(NotImplementedError, lambda: solver._matches())
    raises(NotImplementedError, lambda: solver._get_general_solution())

    # This ODE can not be solved by the FirstLinear solver. Here we test that
    # it does not match and the asking for a general solution gives
    # ODEMatchError

    problem = SingleODEProblem(f(x).diff(x) + f(x)*f(x), f(x), x)

    solver = FirstLinear(problem)
    raises(ODEMatchError, lambda: solver.get_general_solution())

    solver = FirstLinear(problem)
    assert solver.matches() is False

    #These are just test for order of ODE

    problem = SingleODEProblem(f(x).diff(x) + f(x), f(x), x)
    assert problem.order == 1

    problem = SingleODEProblem(f(x).diff(x,4) + f(x).diff(x,2) - f(x).diff(x,3), f(x), x)
    assert problem.order == 4


def test_nth_algebraic():
    eqn = f(x) + f(x)*f(x).diff(x)
    solns = [Eq(f(x), exp(x)),
             Eq(f(x), C1*exp(C2*x))]
    solns_final =  _remove_redundant_solutions(eqn, solns, 2, x)
    assert solns_final == [Eq(f(x), C1*exp(C2*x))]

    _ode_solver_test(_get_examples_ode_sol_nth_algebraic())


def test_factorable():
    _ode_solver_test(_get_examples_ode_sol_factorable())


def test_Riccati_special_minus2():
    _ode_solver_test(_get_examples_ode_sol_riccati())


def test_Bernoulli():
    _ode_solver_test(_get_examples_ode_sol_bernoulli())


def test_1st_linear():
    _ode_solver_test(_get_examples_ode_sol_1st_linear())


def test_almost_linear():
   _ode_solver_test(_get_examples_ode_sol_almost_linear())


def test_Liouville_ODE():
    hint = 'Liouville'
    not_Liouville1 = classify_ode(diff(f(x), x)/x + f(x)*diff(f(x), x, x)/2 -
        diff(f(x), x)**2/2, f(x))
    not_Liouville2 = classify_ode(diff(f(x), x)/x + diff(f(x), x, x)/2 -
        x*diff(f(x), x)**2/2, f(x))
    assert hint not in not_Liouville1
    assert hint not in not_Liouville2
    assert hint + '_Integral' not in not_Liouville1
    assert hint + '_Integral' not in not_Liouville2

    _ode_solver_test(_get_examples_ode_sol_liouville())


def test_nth_order_linear_euler_eq_homogeneous():
    x, t, a, b, c = symbols('x t a b c')
    y = Function('y')
    our_hint = "nth_linear_euler_eq_homogeneous"

    eq = diff(f(t), t, 4)*t**4 - 13*diff(f(t), t, 2)*t**2 + 36*f(t)
    assert our_hint in classify_ode(eq)

    eq = a*y(t) + b*t*diff(y(t), t) + c*t**2*diff(y(t), t, 2)
    assert our_hint in classify_ode(eq)

    _ode_solver_test(_get_examples_ode_sol_euler_homogeneous())


def test_nth_order_linear_euler_eq_nonhomogeneous_undetermined_coefficients():
    x, t = symbols('x t')
    a, b, c, d = symbols('a b c d', integer=True)
    our_hint = "nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients"

    eq = x**4*diff(f(x), x, 4) - 13*x**2*diff(f(x), x, 2) + 36*f(x) + x
    assert our_hint in classify_ode(eq, f(x))

    eq = a*x**2*diff(f(x), x, 2) + b*x*diff(f(x), x) + c*f(x) + d*log(x)
    assert our_hint in classify_ode(eq, f(x))

    _ode_solver_test(_get_examples_ode_sol_euler_undetermined_coeff())


def test_nth_order_linear_euler_eq_nonhomogeneous_variation_of_parameters():
    x, t = symbols('x, t')
    a, b, c, d = symbols('a, b, c, d', integer=True)
    our_hint = "nth_linear_euler_eq_nonhomogeneous_variation_of_parameters"

    eq = Eq(x**2*diff(f(x),x,2) - 8*x*diff(f(x),x) + 12*f(x), x**2)
    assert our_hint in classify_ode(eq, f(x))

    eq = Eq(a*x**3*diff(f(x),x,3) + b*x**2*diff(f(x),x,2) + c*x*diff(f(x),x) + d*f(x), x*log(x))
    assert our_hint in classify_ode(eq, f(x))

    _ode_solver_test(_get_examples_ode_sol_euler_var_para())


def _get_examples_ode_sol_euler_homogeneous():
    return {
            'hint': "nth_linear_euler_eq_homogeneous",
            'func': f(x),
            'examples':{
    'euler_hom_01': {
        'eq': Eq(-3*diff(f(x), x)*x + 2*x**2*diff(f(x), x, x), 0),
        'sol': [Eq(f(x), C1 + C2*x**Rational(5, 2))],
    },

    'euler_hom_02': {
        'eq': Eq(3*f(x) - 5*diff(f(x), x)*x + 2*x**2*diff(f(x), x, x), 0),
        'sol': [Eq(f(x), C1*sqrt(x) + C2*x**3)]
    },

    'euler_hom_03': {
        'eq': Eq(4*f(x) + 5*diff(f(x), x)*x + x**2*diff(f(x), x, x), 0),
        'sol': [Eq(f(x), (C1 + C2*log(x))/x**2)]
    },

    'euler_hom_04': {
        'eq': Eq(6*f(x) - 6*diff(f(x), x)*x + 1*x**2*diff(f(x), x, x) + x**3*diff(f(x), x, x, x), 0),
        'sol': [Eq(f(x), C1/x**2 + C2*x + C3*x**3)]
    },

    'euler_hom_05': {
        'eq': Eq(-125*f(x) + 61*diff(f(x), x)*x - 12*x**2*diff(f(x), x, x) + x**3*diff(f(x), x, x, x), 0),
        'sol': [Eq(f(x), x**5*(C1 + C2*log(x) + C3*log(x)**2))]
    },

    'euler_hom_06': {
        'eq': x**2*diff(f(x), x, 2) + x*diff(f(x), x) - 9*f(x),
        'sol': [Eq(f(x), C1*x**-3 + C2*x**3)]
    },

    'euler_hom_07': {
        'eq': sin(x)*x**2*f(x).diff(x, 2) + sin(x)*x*f(x).diff(x) + sin(x)*f(x),
        'sol': [Eq(f(x), C1*sin(log(x)) + C2*cos(log(x)))],
        'XFAIL': ['2nd_power_series_regular','nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients']
    },
    }
    }


def _get_examples_ode_sol_euler_undetermined_coeff():
    return {
            'hint': "nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients",
            'func': f(x),
            'examples':{
    'euler_undet_01': {
        'eq': Eq(x**2*diff(f(x), x, x) + x*diff(f(x), x), 1),
        'sol': [Eq(f(x), C1 + C2*log(x) + log(x)**2/2)]
    },

    'euler_undet_02': {
        'eq': Eq(x**2*diff(f(x), x, x) - 2*x*diff(f(x), x) + 2*f(x), x**3),
        'sol': [Eq(f(x), x*(C1 + C2*x + Rational(1, 2)*x**2))]
    },

    'euler_undet_03': {
        'eq': Eq(x**2*diff(f(x), x, x) - x*diff(f(x), x) - 3*f(x), log(x)/x),
        'sol': [Eq(f(x), (C1 + C2*x**4 - log(x)**2/8 - log(x)/16)/x)]
    },

    'euler_undet_04': {
        'eq': Eq(x**2*diff(f(x), x, x) + 3*x*diff(f(x), x) - 8*f(x), log(x)**3 - log(x)),
        'sol': [Eq(f(x), C1/x**4 + C2*x**2 - Rational(1,8)*log(x)**3 - Rational(3,32)*log(x)**2 - Rational(1,64)*log(x) - Rational(7, 256))]
    },

    'euler_undet_05': {
        'eq': Eq(x**3*diff(f(x), x, x, x) - 3*x**2*diff(f(x), x, x) + 6*x*diff(f(x), x) - 6*f(x), log(x)),
        'sol': [Eq(f(x), C1*x + C2*x**2 + C3*x**3 - Rational(1, 6)*log(x) - Rational(11, 36))]
    },
    }
    }


def _get_examples_ode_sol_euler_var_para():
    return {
            'hint': "nth_linear_euler_eq_nonhomogeneous_variation_of_parameters",
            'func': f(x),
            'examples':{
    'euler_var_01': {
        'eq': Eq(x**2*Derivative(f(x), x, x) - 2*x*Derivative(f(x), x) + 2*f(x), x**4),
        'sol': [Eq(f(x), x*(C1 + C2*x + x**3/6))]
    },

    'euler_var_02': {
        'eq': Eq(3*x**2*diff(f(x), x, x) + 6*x*diff(f(x), x) - 6*f(x), x**3*exp(x)),
        'sol': [Eq(f(x), C1/x**2 + C2*x + x*exp(x)/3 - 4*exp(x)/3 + 8*exp(x)/(3*x) - 8*exp(x)/(3*x**2))]
    },

    'euler_var_03': {
        'eq': Eq(x**2*Derivative(f(x), x, x) - 2*x*Derivative(f(x), x) + 2*f(x), x**4*exp(x)),
        'sol':  [Eq(f(x), x*(C1 + C2*x + x*exp(x) - 2*exp(x)))]
    },

    'euler_var_04': {
        'eq': x**2*Derivative(f(x), x, x) - 2*x*Derivative(f(x), x) + 2*f(x) - log(x),
        'sol': [Eq(f(x), C1*x + C2*x**2 + log(x)/2 + Rational(3, 4))]
    },

    'euler_var_05': {
        'eq': -exp(x) + (x*Derivative(f(x), (x, 2)) + Derivative(f(x), x))/x,
        'sol': [Eq(f(x), C1 + C2*log(x) + exp(x) - Ei(x))]
    },
    }
    }


def _get_examples_ode_sol_bernoulli():
    # Type: Bernoulli, f'(x) + p(x)*f(x) == q(x)*f(x)**n
    return {
            'hint': "Bernoulli",
            'func': f(x),
            'examples':{
    'bernoulli_01': {
        'eq': Eq(x*f(x).diff(x) + f(x) - f(x)**2, 0),
        'sol': [Eq(f(x), 1/(C1*x + 1))],
        'XFAIL': ['separable_reduced']
    },

    'bernoulli_02': {
        'eq': f(x).diff(x) - y*f(x),
        'sol': [Eq(f(x), C1*exp(x*y))]
    },

    'bernoulli_03': {
        'eq': f(x)*f(x).diff(x) - 1,
        'sol': [Eq(f(x), -sqrt(C1 + 2*x)), Eq(f(x), sqrt(C1 + 2*x))]
    },
    }
    }


def _get_examples_ode_sol_riccati():
    # Type: Riccati special alpha = -2, a*dy/dx + b*y**2 + c*y/x +d/x**2
    return {
            'hint': "Riccati_special_minus2",
            'func': f(x),
            'examples':{
    'riccati_01': {
        'eq': 2*f(x).diff(x) + f(x)**2 - f(x)/x + 3*x**(-2),
        'sol': [Eq(f(x), (-sqrt(3)*tan(C1 + sqrt(3)*log(x)/4) + 3)/(2*x))],
    },
    },
    }


def _get_examples_ode_sol_1st_linear():
    # Type: first order linear form f'(x)+p(x)f(x)=q(x)
    return {
            'hint': "1st_linear",
            'func': f(x),
            'examples':{
    'linear_01': {
        'eq': Eq(f(x).diff(x) + x*f(x), x**2),
        'sol': [Eq(f(x), (C1 + x*exp(x**2/2)- sqrt(2)*sqrt(pi)*erfi(sqrt(2)*x/2)/2)*exp(-x**2/2))],
    },
    },
    }


def _get_examples_ode_sol_factorable():
    """ some hints are marked as xfail for examples because they missed additional algebraic solution
    which could be found by Factorable hint. Fact_01 raise exception for
    nth_linear_constant_coeff_undetermined_coefficients"""

    y = Dummy('y')
    return {
            'hint': "factorable",
            'func': f(x),
            'examples':{
    'fact_01': {
        'eq': f(x) + f(x)*f(x).diff(x),
        'sol': [Eq(f(x), 0), Eq(f(x), C1 - x)],
        'XFAIL': ['separable', '1st_exact', '1st_linear', 'Bernoulli', '1st_homogeneous_coeff_best',
        '1st_homogeneous_coeff_subs_indep_div_dep', '1st_homogeneous_coeff_subs_dep_div_indep',
        'lie_group', 'nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients',
        'nth_linear_constant_coeff_variation_of_parameters',
        'nth_linear_euler_eq_nonhomogeneous_variation_of_parameters',
        'nth_linear_constant_coeff_undetermined_coefficients']
    },

    'fact_02': {
        'eq': f(x)*(f(x).diff(x)+f(x)*x+2),
        'sol': [Eq(f(x), (C1 - sqrt(2)*sqrt(pi)*erfi(sqrt(2)*x/2))*exp(-x**2/2)), Eq(f(x), 0)],
        'XFAIL': ['Bernoulli', '1st_linear', 'lie_group']
    },

    'fact_03': {
        'eq': (f(x).diff(x)+f(x)*x**2)*(f(x).diff(x, 2) + x*f(x)),
        'sol':  [Eq(f(x), C1*airyai(-x) + C2*airybi(-x)),Eq(f(x), C1*exp(-x**3/3))]
    },

    'fact_04': {
        'eq': (f(x).diff(x)+f(x)*x**2)*(f(x).diff(x, 2) + f(x)),
        'sol': [Eq(f(x), C1*exp(-x**3/3)), Eq(f(x), C1*sin(x) + C2*cos(x))]
    },

    'fact_05': {
        'eq': (f(x).diff(x)**2-1)*(f(x).diff(x)**2-4),
        'sol': [Eq(f(x), C1 - x), Eq(f(x), C1 + x), Eq(f(x), C1 + 2*x), Eq(f(x), C1 - 2*x)]
    },

    'fact_06': {
        'eq': (f(x).diff(x, 2)-exp(f(x)))*f(x).diff(x),
        'sol': [Eq(f(x), C1)]
    },

    'fact_07': {
        'eq': (f(x).diff(x)**2-1)*(f(x)*f(x).diff(x)-1),
        'sol': [Eq(f(x), C1 - x), Eq(f(x), -sqrt(C1 + 2*x)),Eq(f(x), sqrt(C1 + 2*x)), Eq(f(x), C1 + x)]
    },

    'fact_08': {
        'eq': Derivative(f(x), x)**4 - 2*Derivative(f(x), x)**2 + 1,
        'sol': [Eq(f(x), C1 - x), Eq(f(x), C1 + x)]
    },

    'fact_09': {
        'eq': f(x)**2*Derivative(f(x), x)**6 - 2*f(x)**2*Derivative(f(x),
         x)**4 + f(x)**2*Derivative(f(x), x)**2 - 2*f(x)*Derivative(f(x),
         x)**5 + 4*f(x)*Derivative(f(x), x)**3 - 2*f(x)*Derivative(f(x),
         x) + Derivative(f(x), x)**4 - 2*Derivative(f(x), x)**2 + 1,
        'sol': [Eq(f(x), C1 - x), Eq(f(x), -sqrt(C1 + 2*x)),
           Eq(f(x), sqrt(C1 + 2*x)), Eq(f(x), C1 + x)]
    },

    'fact_10': {
        'eq': x**4*f(x)**2 + 2*x**4*f(x)*Derivative(f(x), (x, 2)) + x**4*Derivative(f(x),
         (x, 2))**2  + 2*x**3*f(x)*Derivative(f(x), x) + 2*x**3*Derivative(f(x),
         x)*Derivative(f(x), (x, 2)) - 7*x**2*f(x)**2 - 7*x**2*f(x)*Derivative(f(x),
         (x, 2)) + x**2*Derivative(f(x), x)**2 - 7*x*f(x)*Derivative(f(x), x) + 12*f(x)**2,
        'sol': [Eq(f(x), C1*besselj(2, x) + C2*bessely(2, x)), Eq(f(x), C1*besselj(sqrt(3),
           x) + C2*bessely(sqrt(3), x))]
    },

    'fact_11': {
        'eq': (f(x).diff(x, 2)-exp(f(x)))*(f(x).diff(x, 2)+exp(f(x))),
        'sol': [], #currently dsolve doesn't return any solution for this example
        'XFAIL': ['factorable']
    },

    #Below examples were added for the issue: https://github.com/sympy/sympy/issues/15889
    'fact_12': {
        'eq': exp(f(x).diff(x))-f(x)**2,
        'sol': [Eq(NonElementaryIntegral(1/log(y**2), (y, f(x))), C1 + x)],
        'XFAIL': ['lie_group'] #It shows not implemented error for lie_group.
    },

    'fact_13': {
        'eq': f(x).diff(x)**2 - f(x)**3,
        'sol': [Eq(f(x), 4/(C1**2 - 2*C1*x + x**2))],
        'XFAIL': ['lie_group'] #It shows not implemented error for lie_group.
    },

    'fact_14': {
        'eq': f(x).diff(x)**2 - f(x),
        'sol': [Eq(f(x), C1**2/4 - C1*x/2 + x**2/4)]
    },

    'fact_15': {
        'eq': f(x).diff(x)**2 - f(x)**2,
        'sol': [Eq(f(x), C1*exp(x)), Eq(f(x), C1*exp(-x))]
    },

    'fact_16': {
        'eq': f(x).diff(x)**2 - f(x)**3,
        'sol': [Eq(f(x), 4/(C1**2 - 2*C1*x + x**2))]
    },
    }
    }



def _get_examples_ode_sol_almost_linear():
    from sympy import Ei
    A = Symbol('A', positive=True)
    f = Function('f')
    d = f(x).diff(x)

    return {
            'hint': "almost_linear",
            'func': f(x),
            'examples':{
    'almost_lin_01': {
        'eq': x**2*f(x)**2*d + f(x)**3 + 1,
        'sol': [Eq(f(x), (C1*exp(3/x) - 1)**Rational(1, 3)),
        Eq(f(x), (-1 - sqrt(3)*I)*(C1*exp(3/x) - 1)**Rational(1, 3)/2),
        Eq(f(x), (-1 + sqrt(3)*I)*(C1*exp(3/x) - 1)**Rational(1, 3)/2)],

    },

    'almost_lin_02': {
        'eq': x*f(x)*d + 2*x*f(x)**2 + 1,
        'sol': [Eq(f(x), -sqrt((C1 - 2*Ei(4*x))*exp(-4*x))), Eq(f(x), sqrt((C1 - 2*Ei(4*x))*exp(-4*x)))]
    },

    'almost_lin_03': {
        'eq':  x*d + x*f(x) + 1,
        'sol': [Eq(f(x), (C1 - Ei(x))*exp(-x))]
    },

    'almost_lin_04': {
        'eq': x*exp(f(x))*d + exp(f(x)) + 3*x,
        'sol': [Eq(f(x), log(C1/x - x*Rational(3, 2)))],
    },

    'almost_lin_05': {
        'eq': x + A*(x + diff(f(x), x) + f(x)) + diff(f(x), x) + f(x) + 2,
        'sol': [Eq(f(x), (C1 + Piecewise(
        (x, Eq(A + 1, 0)), ((-A*x + A - x - 1)*exp(x)/(A + 1), True)))*exp(-x))],
    },
    }
    }


def _get_examples_ode_sol_liouville():
    return {
            'hint': "Liouville",
            'func': f(x),
            'examples':{
    'liouville_01': {
        'eq': diff(f(x), x)/x + diff(f(x), x, x)/2 - diff(f(x), x)**2/2,
        'sol': [Eq(f(x), log(x/(C1 + C2*x)))],

    },

    'liouville_02': {
        'eq': diff(x*exp(-f(x)), x, x),
        'sol': [Eq(f(x), log(x/(C1 + C2*x)))]
    },

    'liouville_03': {
        'eq':  ((diff(f(x), x)/x + diff(f(x), x, x)/2 - diff(f(x), x)**2/2)*exp(-f(x))/exp(f(x))).expand(),
        'sol': [Eq(f(x), log(x/(C1 + C2*x)))]
    },

    'liouville_04': {
        'eq': diff(f(x), x, x) + 1/f(x)*(diff(f(x), x))**2 + 1/x*diff(f(x), x),
        'sol': [Eq(f(x), -sqrt(C1 + C2*log(x))), Eq(f(x), sqrt(C1 + C2*log(x)))],
    },

    'liouville_05': {
        'eq': x*diff(f(x), x, x) + x/f(x)*diff(f(x), x)**2 + x*diff(f(x), x),
        'sol': [Eq(f(x), -sqrt(C1 + C2*exp(-x))), Eq(f(x), sqrt(C1 + C2*exp(-x)))],
    },

    'liouville_06': {
        'eq': Eq((x*exp(f(x))).diff(x, x), 0),
        'sol': [Eq(f(x), log(C1 + C2/x))],
    },
    }
    }


def _get_examples_ode_sol_nth_algebraic():
    M, m, r, t = symbols('M m r t')
    phi = Function('phi')
    # This one needs a substitution f' = g.
    # 'algeb_12': {
    #     'eq': -exp(x) + (x*Derivative(f(x), (x, 2)) + Derivative(f(x), x))/x,
    #     'sol': [Eq(f(x), C1 + C2*log(x) + exp(x) - Ei(x))],
    # },
    return {
            'hint': "nth_algebraic",
            'func': f(x),
            'examples':{
    'algeb_01': {
        'eq': f(x) * f(x).diff(x) * f(x).diff(x, x) * (f(x) - 1) * (f(x).diff(x) - x),
        'sol': [Eq(f(x), C1 + x**2/2), Eq(f(x), C1 + C2*x)]
    },

    'algeb_02': {
        'eq': f(x) * f(x).diff(x) * f(x).diff(x, x) * (f(x) - 1),
        'sol': [Eq(f(x), C1 + C2*x)]
    },

    'algeb_03': {
        'eq': f(x) * f(x).diff(x) * f(x).diff(x, x),
        'sol': [Eq(f(x), C1 + C2*x)]
    },

    'algeb_04': {
        'eq': Eq(-M * phi(t).diff(t),
         Rational(3, 2) * m * r**2 * phi(t).diff(t) * phi(t).diff(t,t)),
        'sol': [Eq(phi(t), C1), Eq(phi(t), C1 + C2*t - M*t**2/(3*m*r**2))],
        'func': phi(t)
    },

    'algeb_05': {
        'eq': (1 - sin(f(x))) * f(x).diff(x),
        'sol': [Eq(f(x), C1)],
        'XFAIL': ['separable']  #It raised exception.
    },

    'algeb_06': {
        'eq': (diff(f(x)) - x)*(diff(f(x)) + x),
        'sol': [Eq(f(x), C1 - x**2/2), Eq(f(x), C1 + x**2/2)]
    },

    'algeb_07': {
        'eq': Eq(Derivative(f(x), x), Derivative(g(x), x)),
        'sol': [Eq(f(x), C1 + g(x))],
    },

    'algeb_08': {
        'eq': f(x).diff(x) - C1,   #this example is from issue 15999
        'sol': [Eq(f(x), C1*x + C2)],
    },

    'algeb_09': {
        'eq': f(x)*f(x).diff(x),
        'sol': [Eq(f(x), C1)],
    },

    'algeb_10': {
        'eq': (diff(f(x)) - x)*(diff(f(x)) + x),
        'sol': [Eq(f(x), C1 - x**2/2), Eq(f(x), C1 + x**2/2)],
    },

    'algeb_11': {
        'eq': f(x) + f(x)*f(x).diff(x),
        'sol': [Eq(f(x), 0), Eq(f(x), C1 - x)],
        'XFAIL': ['separable', '1st_exact', '1st_linear', 'Bernoulli', '1st_homogeneous_coeff_best',
         '1st_homogeneous_coeff_subs_indep_div_dep', '1st_homogeneous_coeff_subs_dep_div_indep',
         'lie_group', 'nth_linear_constant_coeff_undetermined_coefficients',
         'nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients',
         'nth_linear_constant_coeff_variation_of_parameters',
         'nth_linear_euler_eq_nonhomogeneous_variation_of_parameters']
         #nth_linear_constant_coeff_undetermined_coefficients raises exception rest all of them misses a solution.
    },

    'algeb_12': {
        'eq': Derivative(x*f(x), x, x, x),
        'sol': [Eq(f(x), (C1 + C2*x + C3*x**2) / x)],
        'XFAIL': ['nth_algebraic']  # It passes only when prep=False is set in dsolve.
    },

    'algeb_13': {
        'eq': Eq(Derivative(x*Derivative(f(x), x), x)/x, exp(x)),
        'sol': [Eq(f(x), C1 + C2*log(x) + exp(x) - Ei(x))],
        'XFAIL': ['nth_algebraic']  # It passes only when prep=False is set in dsolve.
    },
    }
    }


def _get_all_examples():
    all_solvers = [_get_examples_ode_sol_euler_homogeneous(),
    _get_examples_ode_sol_euler_undetermined_coeff(),
    _get_examples_ode_sol_euler_var_para(),
    _get_examples_ode_sol_factorable(),
    _get_examples_ode_sol_bernoulli(),
    _get_examples_ode_sol_nth_algebraic(),
    _get_examples_ode_sol_riccati(),
    _get_examples_ode_sol_1st_linear(),
    _get_examples_ode_sol_almost_linear(),
    ]

    all_examples = []
    for solver in all_solvers:
        for example in solver['examples']:
            temp = {
                'hint': solver['hint'],
                'func': solver['examples'][example].get('func',solver['func']),
                'eq': solver['examples'][example]['eq'],
                'sol': solver['examples'][example]['sol'],
                'XFAIL': solver['examples'][example].get('XFAIL',[]),
                'example_name': example,
            }
            all_examples.append(temp)
    return all_examples
