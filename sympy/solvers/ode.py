"""
This module contains dsolve() and different helper functions that it uses.

dsolve() solves ordinary differential equations. See the docstring on
the various functions for their uses. Note that partial differential
equations support is in pde.py.  Note that ode_hint() functions have
docstrings describing their various methods, but they are intended
for internal use.  Use dsolve(ode, func, hint=hint) to solve an ode
using a specific hint.  See also the docstring on dsolve().

== Functions in this module ==

    These are the user functions in this module:

    - dsolve() - Solves ODEs.
    - classify_ode() - Classifies ODEs into possible hints for dsolve().
    - checkodesol() - Checks if an equation is the solution to an ODE.
    - ode_order() - Returns the order (degree) of an ODE.
    - homogeneous_order() - Returns the homogeneous order of an
      expression.

    See also the docstrings of these functions.  There are also quite
    a few functions that are not imported into the global namespace.
    See the docstrings of those functions for more info.

== Solving methods currently implemented ==

The following methods are implemented for solving ordinary differential
equations.  See the docstrings of the various ode_hint() functions for
more information on each (run help(ode)):
    - 1st order separable differential equations
    - 1st order differential equations whose coefficients or dx and dy
      are functions homogeneous of the same order.
    - 1st order exact differential equations.
    - 1st order linear differential equations
    - 1st order Bernoulli differential equations.
    - 2nd order Liouville differential equations.
    - nth order linear homogeneous differential equation with constant
      coefficients.
    - nth order linear inhomogeneous differential equation with constant
      coefficients using the method of undetermined coefficients.
    - nth order linear inhomogeneous differential equation with constant
      coefficients using the method of variation of parameters.

== Philosophy behind this module ==

This module is designed to make it easy to add new ODE solving methods
without having to mess with the solving code for other methods.  The
idea is that there is a classify_ode() function, which takes in an ODE
and tells you what hints, if any, will solve the ODE.  It does this
without attempting to solve the ODE, so it is fast.  Each solving method
is a hint, and it has its own function, named ode_hint.  That function
takes in the ODE and any match expression gathered by classify_ode and
returns a solved result.  If this result has any integrals in it, the
ode_hint function will return an unevaluated Integral class. dsolve(),
which is the user wrapper function around all of this, will then call
odesimp() on the result, which, among other things, will attempt to
solve the equation for the dependent variable (the function we are
solving for), simplify the arbitrary constants in the expression, and
evaluate any integrals, if the hint allows it.

== How to add new solution methods ==

If you have an ODE that you want dsolve() to be able to solve, try to
avoid adding special case code here.  Instead, try finding a general
method that will solve your ODE, as well as others.  This way, the ode
module will become more robust, and unhindered by special case hacks.
WolphramAlpha and Maple's DETools[odeadvisor] function are two resources
you can use to classify a specific ODE.  It is also better for a method
to work with an nth order ODE instead of only with specific orders, if
possible.

To add a new method, there are a few things that you need to do.  First,
you need a hint name for your method.  Try to name your hint so that it
is unambiguous with all other methods, include ones that may not be
implemented yet.  If your method uses integrals, also include a
"hint_Integral" hint.  If there is more than one way to solve ODEs with
your method, include a hint for each one, as well as a "hint_best" hint.
Your ode_hint_best() function should choose the best using
compare_ode_sol.  See ode_1st_homogeneous_coeff_best(), for example. The
function that uses your method will be called ode_hint(), so the hint
must only use characters that are allowed in a Python function name
(alphanumeric characters and the underscore '_' character).  Include a
function for every hint, except for "_Integral" hints (dsolve() takes
care of those automatically).  Hint names should be all lowercase,
unless a word is commonly capitalized (such as Integral or Bernoulli).
If you have a hint that you do not want to run with "all_Integral" that
doesn't have an "_Integral" counterpart (such as a best hint that would
defeat the purpose of "all_Integral"), you will need to remove it
manually in the dsolve() code.  See also the classify_ode() docstring
for guidelines on writing a hint name.

Determine **in general** how the solutions returned by your method
compare with other methods that can potentially solve the same ODEs.
Then, put your hints in the allhints tuple in the order that they should
be called.  The ordering of this tuple determines which hints are
default. Note that exceptions are ok, because it is easy for the user to
choose individual hints with dsolve().  In general, "_Integral" variants
should go at the end of the list, and "_best" variants should go before
the various hints they apply to.  For example, the
"undetermined_coefficients" hint comes before the
"variation_of_parameters" hint because, even though variation of
parameters can be used to solve all ODEs that undetermined coefficients
can solve and more, undetermined coefficients generally returns cleaner
results for the ODEs that it can solve than variation of parameters
does, and it does not require integration, so it is much faster.

Next, you need to have a match expression or a function that matches the
type of the ODE, which you should put in classify_ode() (if the match
function is more than just a few lines, like
_undetermined_coefficients_match(), it should go outside of
classify_ode()).  It should match the ODE without solving for it as much
as possible, so that classify_ode() remains fast and is not hindered by
bugs in solving code.  Be sure to consider corner cases. For example, if
your solution method involves dividing by something, make sure you
exclude the case where that division will be 0.

In most cases, the matching of the ODE will also give you the various
parts that you need to solve it. You should put that in a dictionary
(.match() will do this for you), and add that as matching_hints['hint']
= matchdict in the relevant part of classify_ode.  classify_ode will
then send this to dsolve(), which will send it to your function as the
match argument. Your function should be named ode_hint(eq, func, order,
match). If you need to send more information, put it in the dictionary
match.  For example, if you used a dummy variable in classify_ode to
match your expression, you will ned to pass it to your function using
the match dict to access it.  You can access the independent variable
using func.args[0], and the dependent variable (function to solve for)
as func.func.  If, while trying to solve the ODE, you find that you
cannot, raise NotImplementedError.  dsolve() will catch this error with
the "all" meta-hint, rather than causing the whole routine to fail.

Add a docstring to your function that describes the method employed.
Like with anything else in SymPy, you will need to add a doctest to the
docstring, in addition to real tests in test_ode.py.  Try to maintain
consistency with the other hint functions' docstrings.  Add your method
to the list at the top of this docstring.  Also, add your method to
ode.txt in the docs/src directory, so that the Sphinx docs will pull its
docstring into the main SymPy documentation.

If your solution method involves integrating, use C.Integral() instead
of integrate().  This allows the user to bypass hard/slow integration by
using the "_Integral" variant of your hint.  In most cases, calling
.doit() will integrate your solution.  If this is not the case, you will
need to write special code in _handle_Integral().  Arbitrary constants
should be symbols named C1, C2, and so on.  All solution methods should
return an equality instance.  If you need an arbitrary number of
arbitrary constants, you can use constants =
numbered_symbols(prefix='C', function=Symbol, start=1).  If it is
possible to solve for the dependent function in a general way, do so.
Otherwise, do as best as you can, but do not call solve in your
ode_hint() function.  odesimp() will attempt to solve the solution for
you, so you do not need to do that. Lastly, if your ODE has a common
simplification that can be applied to your solutions, you can add a
special case in odesimp() for it.  For example, solutions returned from
the "1st_homogeneous_coeff" hints often have many log() terms, so
odesimp() calls logcombine() on them (it also helps to write the
arbitrary constant as log(C1) instead of C1 in this case).  Also
consider common ways that you can rearrange your solution to have
constantsimp() take better advantage of it.  It is better to put
simplification in odesimp() than in your method, because it can then be
turned off with the simplify flag in dsolve(). If you have any
extraneous simplification in your function, be sure to only run it using
"if match.get('simplify', True):", especially if it can be slow or if it
can reduce the domain of the solution.

Feel free to refactor existing hints to avoid duplicating code or
creating inconsistencies.  If you can show that your method exactly
duplicates an existing method, including in the simplicity and speed of
obtaining the solutions, then you can remove the old, less general
method.  The existing code is tested extensively in test_ode.py, so if
anything is broken, one of those tests will surely fail.

"""
from sympy.core.basic import Add, Basic, C, Mul, Pow
from sympy.core.function import Derivative, diff, expand_mul, Function
from sympy.core.multidimensional import vectorize
from sympy.core.relational import Equality, Eq
from sympy.core.symbol import Symbol, Wild
from sympy.core.sympify import sympify

from sympy.functions import cos, exp, im, log, re, sin
from sympy.matrices import wronskian
from sympy.polys import RootsOf, discriminant, RootOf
from sympy.simplify import collect, logcombine, powsimp, separatevars, \
    simplify, trigsimp
from sympy.solvers import solve

from sympy.utilities import numbered_symbols

import sympy.solvers
# This is a list of hints in the order that they should be applied.  That means
# that, in general, hints earlier in the list should produce simpler results
# than those later for ODEs that fit both.  This is just based on my own
# empirical observations, so if you find that *in general*, a hint later in
# the list is better than one before it, fell free to modify the list.  Note
# however that you can easily override the hint used in dsolve() for a specific ODE
# (see the docstring).  In general, "_Integral" hints should be grouped
# at the end of the list, unless there is a method that returns an unevaluatable
# integral most of the time (which should surely go near the end of the list
# anyway).
# "default", "all", "best", and "all_Integral" meta-hints should not be
# included in this list, but "_best" and "_Integral" hints should be included.
allhints = ("separable", "1st_exact", "1st_linear", "Bernoulli",
"1st_homogeneous_coeff_best", "1st_homogeneous_coeff_subs_indep_div_dep",
"1st_homogeneous_coeff_subs_dep_div_indep", "nth_linear_constant_coeff_homogeneous",
"nth_linear_constant_coeff_undetermined_coefficients",
"nth_linear_constant_coeff_variation_of_parameters",
"Liouville", "separable_Integral", "1st_exact_Integral", "1st_linear_Integral",
"Bernoulli_Integral", "1st_homogeneous_coeff_subs_indep_div_dep_Integral",
"1st_homogeneous_coeff_subs_dep_div_indep_Integral",
"nth_linear_constant_coeff_variation_of_parameters_Integral",
"Liouville_Integral")


def dsolve(eq, func, hint="default", simplify=True, **kwargs):
    """
    Solves any (supported) kind of ordinary differential equation.

    == Usage ==

        dsolve(eq, f(x), hint) -> Solve ordinary differential equation
        eq for function f(x), using method hint.


    == Details ==

        eq can be any supported ordinary differential equation (see the
            ode docstring for supported methods).  This can either be an
            Equality, or an expression, which is assumed to be equal to 0.

        f(x) is a function of one variable whose derivatives in that
            variable make up the ordinary differential equation eq.

        hint is the solving method that you want dsolve to use.  Use
            classify_ode(eq, f(x)) to get all of the possible hints for
            an ODE.  The default hint, 'default', will use whatever
            hint is returned first by classify_ode().  See Hints below
            for more options that you can use for hint.
        simplify enables simplification by odesimp().  See its docstring
            for more information.  Turn this off, for example, to
            disable solving of solutions for func or simplification of
            arbitrary constants.  It will still integrate with this hint.
            Note that the solution may contain more arbitrary constants
            than the order of the ODE with this option enabled.

    == Hints ==

        Aside from the various solving methods, there are also some
        meta-hints that you can
        pass to dsolve():

        "default":
                This uses whatever hint is returned first by
                classify_ode(). This is the default argument to dsolve().

        "all":
                To make dsolve apply all relevant classification hints,
                use dsolve(ODE, func, hint="all").  This will return a
                dictionary of hint:solution terms.  If a hint causes
                dsolve to raise NotImplementedError, value of that
                hint's key will be the exception object raised.  The
                dictionary will also include some special keys:

                - order: The order of the ODE.  See also ode_order().
                - best: The simplest hint; what would be returned by
                  "best" below.
                - best_hint: The hint that would produce the solution
                  given by 'best'.  If more than one hint produces the
                  best solution, the first one in the tuple returned by
                  classify_ode() is chosen.
                - default: The solution that would be returned by
                  default.  This is the one produced by the hint that
                  appears first in the tuple returned by classify_ode().

        "all_Integral":
                This is the same as "all", except if a hint also has a
                corresponding "_Integral" hint, it only returns the
                "_Integral" hint.  This is useful if "all" causes
                dsolve() to hang because of a difficult or impossible
                integral.  This meta-hint will also be much faster than
                "all", because integrate() is an expensive routine.

        "best":
                To have dsolve() try all methods and return the simplest
                one.  This takes into account whether the solution is
                solvable in the function, whether it contains any
                Integral classes (i.e. unevaluatable integrals), and
                which one is the shortest in size.

        See also the classify_ode() docstring for more info on hints,
        and the ode docstring for a list of all supported hints.


    == Tips ==
        - You can declare the derivative of an unknown function this way:
            >>> from sympy import *
            >>> x = Symbol('x') # x is the independent variable
            >>> f = Function("f")(x) # f is a function of x
            >>> # f_ will be the derivative of f with respect to x
            >>> f_ = Derivative(f, x)

        - See test_ode.py for many tests, which serves also as a set of
          examples for how to use dsolve().
        - dsolve always returns an Equality class.  If possible, it
          solves the solution explicitly for the function being solved
          for. Otherwise, it returns an implicit solution.
        - Arbitrary constants are symbols named C1, C2, and so on.
        - Because all solutions should be mathematically equivalent,
          some hints may return the exact same result for an ODE.  Often,
          though, two different hints will return the same solution
          formatted differently.  The two should be equivalent.  Also
          note that sometimes the values of the arbitrary constants in
          two different solutions may not be the same, because one
          constant may have "absorbed" other constants into it.
        - Do help(ode.ode_hintname) to get help more information on a
          specific hint, where hintname is the name of a hint without
          "_Integral".

    == Examples ==

        >>> from sympy import *
        >>> x = Symbol('x')
        >>> f = Function('f')
        >>> dsolve(Derivative(f(x),x,x)+9*f(x), f(x))
        f(x) == C1*sin(3*x) + C2*cos(3*x)
        >>> dsolve(sin(x)*cos(f(x)) + cos(x)*sin(f(x))*f(x).diff(x), f(x),
        ...     hint='separable')
        -log(1 - sin(f(x))**2)/2 == C1 + log(1 - sin(x)**2)/2
        >>> dsolve(sin(x)*cos(f(x)) + cos(x)*sin(f(x))*f(x).diff(x), f(x),
        ...     hint='1st_exact')
        f(x) == acos(C1/cos(x))
        >>> dsolve(sin(x)*cos(f(x)) + cos(x)*sin(f(x))*f(x).diff(x), f(x),
        ... hint='best')
        f(x) == acos(C1/cos(x))
        >>> # Note that even though separable is the default, 1st_exact produces
        >>> # a simpler result in this case.

    """
    # TODO: Implement initial conditions
    if isinstance(eq, Equality):
        if eq.rhs != 0:
            return dsolve(eq.lhs-eq.rhs, func, hint=hint, simplify=simplify, **kwargs)
        eq = eq.lhs

    # Magic that should only be used internally.  Prevents classify_ode from
    # being called more than it needs to be by passing its results through
    # recursive calls.
    if not kwargs.has_key('classify') or kwargs['classify']:
        hints = classify_ode(eq, func, dict=True)
    else:
        if kwargs.has_key('hints'):
            hints = kwargs['hints']
        else:
            hints = {'default': hint, hint: kwargs['match'], 'order': kwargs['order']}

    if not hints['default']:
        # classify_ode will set hints['default'] to None if no hints match.
        raise NotImplementedError("dsolve: Cannot solve " + str(eq))

    if hints['order'] == 0:
        raise ValueError(str(eq) + "is not a differential equation in " + str(func))

    if hint == 'default':
        return dsolve(eq, func, hint=hints['default'], simplify=simplify, classify=False,
        order=hints['order'], match=hints[hints['default']])
    elif hint == 'best':
        return dsolve(eq, func, hint='all', simplify=simplify, classify=False,
            order=hints['order'], hints = hints)['best']
    elif hint in ('all', 'all_Integral'):
        retdict = {}
        failedhints = {}
        gethints = set(hints) - set(['order', 'default', 'ordered_hints'])
        if hint == 'all_Integral':
            for i in hints:
                if i[-9:] == '_Integral':
                    gethints.remove(i[:-9])
            # special case
            if "1st_homogeneous_coeff_best" in gethints:
                gethints.remove("1st_homogeneous_coeff_best")
        for i in gethints:
            try:
                sol = dsolve(eq, func, hint=i, simplify=simplify, classify=False,
                   order=hints['order'], match=hints[i])
            except NotImplementedError, detail: # except NotImplementedError as detail:
                failedhints[i] = detail
            else:
                retdict[i] = sol
        retdict['best'] = sorted(retdict.values(), cmp=lambda x, y:\
            compare_ode_sol(x, y, func))[0]
        for i in hints['ordered_hints']:
            if retdict['best'] == retdict.get(i, None):
                retdict['best_hint'] = i
                break
        retdict['default'] = hints['default']
        retdict['order'] = sympify(hints['order'])
        retdict.update(failedhints)
        return retdict
    elif hint not in allhints:# and hint not in ('default', 'ordered_hints'):
        raise ValueError("Hint not recognized: " + hint)
    elif hint not in hints:
        raise ValueError("ODE " + str(eq) + " does not match hint " + hint)
    elif hint[-9:] == '_Integral':
        solvefunc = globals()['ode_' + hint[:-9]]
    else:
        solvefunc = globals()['ode_' + hint] # convert the string into a function
    # odesimp() will attempt to integrate, if necessary, apply constantsimp(),
    # attempt to solve for func, and apply any other hint specific simplifications
    if simplify:
        return odesimp(solvefunc(eq, func, order=hints['order'],
            match=hints[hint]), func, hints['order'], hint)
    else:
        # We still want to integrate (you can disable it separately with the hint)
        r = hints[hint]
        r['simplify'] = False # Some hints can take advantage of this option
        return _handle_Integral(solvefunc(eq, func, order=hints['order'],
            match=hints[hint]), func, hints['order'], hint)


def classify_ode(eq, func, dict=False):
    """
    Returns a tuple of possible dsolve() classifications for an ODE.

    The tuple is ordered so that first item is the classification that
    dsolve() uses to solve the ODE by default.  In general,
    classifications at the near the beginning of the list will produce
    better solutions faster than those near the end, thought there are
    always exceptions.  To make dsolve use a different classification,
    use dsolve(ODE, func, hint=<classification>).  See also the dsolve()
    docstring for different meta-hints you can use.

    If dict is true, classify_ode() will return a dictionary of
    hint:match expression terms. This is indendened for internal use by
    dsolve().  Note that because dictionaries are ordered arbitrarily,
    this will most likely not be in the same order as the tuple.

    You can get help on different hints by doing help(ode.ode_hintname),
    where hintname is the name of the hint without "_Integral".

    == Notes on Hint Names ==

    === "_Integral" ===

        If a classification has "_Integral" at the end, it will return
        the expression with an unevaluated Integral class in it.  Note
        that a hint may do this anyway if integrate() cannot do the
        integral, though just using an "_Integral" will do so much
        faster.  Indeed, an "_Integral" hint will always be faster than
        its corresponding hint without "_Integral" because integrate()
        is an expensive routine.  If dsolve() hangs, it is probably
        because integrate() is hanging on a tough or impossible
        integral.  Try using an "_Integral" hint or "all_Integral" to
        get it return something.

        Note that some hints do not have "_Integral" counterparts.  This
        is because integrate() is not used in solving the ODE for those
        method. For example, nth order linear homogeneous ODEs with
        constant coefficients do not require integration to solve, so
        there is no "nth_linear_homogeneous_constant_coeff_Integrate"
        hint. You can easliy evaluate any unevaluated Integrals in an
        expression by doing expr.doit().

    === Ordinals ===

        Some hints contain an ordinal such as "1st_linear".  This is to
        help differentiate them from other hints, as well as from other
        methods that may not be implemented yet. If a hint has "nth" in
        it, such as the "nth_linear" hints, this means that the method
        used to applies to ODEs of any order.

    === "indep" and "dep" ===

        Some hints contain the words "indep" or "dep".  These reference
        the independent variable and the dependent function,
        respectively. For example, if an ODE is in terms of f(x), then
        "indep" will refer to x and "dep" will refer to f.

    === "subs" ===

        If a hints has the word "subs" in it, it means the the ODE is
        solved by substituting the expression given after the word
        "subs" for a single dummy variable.  This is usually in terms of
        "indep" and "dep" as above.  The substituted expression will be
        written only in characters allowed for names of Python objects,
        meaning operators will be spelled out.  For example, indep/dep
        will be written as indep_div_dep.

    === "coeff" ===

        The word "coeff" in a hint refers to the coefficients of
        something in the ODE, usually of the derivative terms.  See the
        docstring for the individual methods for more info (help(ode)).
        This is contrast to "coefficients", as in
        "undetermined_coefficients", which refers to the common name of
        a method.

    === "_best" ===

        Methods that have more than one fundamental way to solve will
        have a hint for each sub-method and a "_best"
        meta-classification. This will evaluate all hints and return the
        best, using the same considerations as the normal "best"
        meta-hint.


    == Examples ==
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> f = Function('f')
        >>> classify_ode(Eq(f(x).diff(x), 0), f(x))
        ('separable', '1st_linear', '1st_homogeneous_coeff_best',
        '1st_homogeneous_coeff_subs_indep_div_dep',
        '1st_homogeneous_coeff_subs_dep_div_indep',
        'nth_linear_constant_coeff_homogeneous', 'separable_Integral',
        '1st_linear_Integral',
        '1st_homogeneous_coeff_subs_indep_div_dep_Integral',
        '1st_homogeneous_coeff_subs_dep_div_indep_Integral')
        >>> classify_ode(f(x).diff(x, 2) + 3*f(x).diff(x) + 2*f(x) - 4, f(x))
        ('nth_linear_constant_coeff_undetermined_coefficients',
        'nth_linear_constant_coeff_variation_of_parameters',
        'nth_linear_constant_coeff_variation_of_parameters_Integral')

    """
    from sympy.core import S
    from sympy.utilities import all
    if len(func.args) != 1:
        raise ValueError("dsolve() and classify_ode() only work with functions " + \
            "of one variable")
    x = func.args[0]
    f = func.func
    y = Symbol('y', dummy=True)
    if isinstance(eq, Equality):
        if eq.rhs != 0:
            return classify_ode(eq.lhs-eq.rhs, func)
        eq = eq.lhs
    # Collect diff(f(x),x) terms so that match will work correctly
    # collect() needs to be improved, as this doesn't always work.
    eq = collect(eq, f(x).diff(x))
    order = ode_order(eq, f(x))

    # hint:matchdict or hint:(tuple of matchdicts)
    # Also will contain "default":<default hint> and "order":order items.
    matching_hints = {"order": order}

    a = Wild('a', exclude=[f(x)])
    b = Wild('b', exclude=[f(x)])
    c = Wild('c', exclude=[f(x)])
    d = Wild('d', exclude=[f(x).diff(x), f(x).diff(x, 2)])
    e = Wild('e', exclude=[f(x).diff(x)])
    k = Wild('k', exclude=[f(x).diff(x)])
    n = Wild('n', exclude=[f(x)])

    if order == 1:
        # We can save a lot of time by skipping these if the ODE isn't 1st order

        # Linear case: a(x)*y'+b(x)*y+c(x) == 0
        r = eq.match(a*diff(f(x),x) + b*f(x) + c)
        if r:
            r['a'] = a
            r['b'] = b
            r['c'] = c
            matching_hints["1st_linear"] = r
            matching_hints["1st_linear_Integral"] = r

        # Bernoulli case: a(x)*y'+b(x)*y+c(x)*y**n == 0
        r = eq.match(a*diff(f(x),x) + b*f(x) + c*f(x)**n)
        if r and r[c] != 0 and r[n] != 1: # See issue 1577
            r['a'] = a
            r['b'] = b
            r['c'] = c
            r['n'] = n
            matching_hints["Bernoulli"] = r
            matching_hints["Bernoulli_Integral"] = r

        # This match is used for several cases below.
        r = eq.match(d+e*diff(f(x),x))
        if r:
            r['d'] = d
            r['e'] = e
            r['y'] = y
            r[d] = r[d].subs(f(x),y)
            r[e] = r[e].subs(f(x),y)

            # Separable Case: y' == P(y)*Q(x)
            r[d] = separatevars(r[d])
            r[e] = separatevars(r[e])
            # m1[coeff]*m1[x]*m1[y] + m2[coeff]*m2[x]*m2[y]*y'
            m1 = separatevars(r[d], dict=True, symbols=(x, y))
            m2 = separatevars(r[e], dict=True, symbols=(x, y))
            if m1 and m2:
                r1 = {'m1':m1, 'm2':m2, 'y':y}
                matching_hints["separable"] = r1
                matching_hints["separable_Integral"] = r1

            # Exact Differential Equation: P(x,y)+Q(x,y)*y'=0 where dP/dy == dQ/dx
            if simplify(r[d].diff(y)) == simplify(r[e].diff(x)) and r[d] != 0:
                matching_hints["1st_exact"] = r
                matching_hints["1st_exact_Integral"] = r

            # First order equation with homogeneous coefficients:
            # dy/dx == F(y/x) or dy/dx == F(x/y)
            ordera = homogeneous_order(r[d], x, y)
            orderb = homogeneous_order(r[e], x, y)
            if ordera == orderb and ordera != None:
                # u1=y/x and u2=x/y
                u1 = Symbol('u1', dummy=True)
                u2 = Symbol('u2', dummy=True)
                if simplify((r[d]+u1*r[e]).subs({x:1, y:u1})) != 0:
                    matching_hints["1st_homogeneous_coeff_subs_dep_div_indep"] = r
                    matching_hints["1st_homogeneous_coeff_subs_dep_div_indep_Integral"] = r
                if simplify((r[e]+u2*r[d]).subs({x:u2, y:1})) != 0:
                    matching_hints["1st_homogeneous_coeff_subs_indep_div_dep"] = r
                    matching_hints["1st_homogeneous_coeff_subs_indep_div_dep_Integral"] = r
                if matching_hints.has_key("1st_homogeneous_coeff_subs_dep_div_indep") \
                and matching_hints.has_key("1st_homogeneous_coeff_subs_indep_div_dep"):
                    matching_hints["1st_homogeneous_coeff_best"] = r

    if order == 2:
        # Liouville ODE f(x).diff(x, 2) + g(f(x))*(f(x).diff(x, 2))**2 + h(x)*f(x).diff(x)
        # See Goldstein and Braun, "Advanced Methods for the Solution of
        # Differential Equations", pg. 98
        s = d*f(x).diff(x, 2) + e*f(x).diff(x)**2 + k*f(x).diff(x)
        r = eq.match(s)
        if r and r[d] != 0:
            y = Symbol('y', dummy=True)
            g = simplify(r[e]/r[d]).subs(f(x), y)
            h = simplify(r[k]/r[d])
            if h.has(f(x)) or g.has(x):
                pass
            else:
                r = {'g':g, 'h':h, 'y':y}
                matching_hints["Liouville"] = r
                matching_hints["Liouville_Integral"] = r


    if order > 0:
        # nth order linear ODE

        # I used to use the below match, but bugs in match would prevent
        # it from matching in all cases, so I wrote _match_nth_linear instead.
        # See issues 1429 and 1601.
        # a_n(x)y^(n) + ... + a_1(x)y' + a_0(x)y = F(x)
#        j = 0
#        s = S(0)
#        wilds = []
#        # Build a match expression for a nth order linear ode
#        for i in numbered_symbols(prefix='a', function=Wild, exclude=[f(x)]):
#            if j == order+1:
#                break
#            wilds.append(i)
#            s += i*f(x).diff(x,j)
#            j += 1
#        s += b

#        r = eq.match(s)

        r = _nth_linear_match(eq, func, order) # Alternate matching function

        # Constant coefficient case (a_i is constant for all i)
        if r and all([not r[i].has(x) for i in range(order + 1)]):
            # Inhomogeneous case: F(x) is not identically 0
            if r['b']:
                undetcoeff = _undetermined_coefficients_match(r['b'], x)
                matching_hints["nth_linear_constant_coeff_variation_of_parameters"] = r
                matching_hints["nth_linear_constant_coeff_variation_of_parameters" + \
                    "_Integral"] = r
                if undetcoeff['test']:
                    r['trialset'] = undetcoeff['trialset']
                    matching_hints["nth_linear_constant_coeff_undetermined_" + \
                        "coefficients"] = r
            # Homogeneous case: F(x) is identically 0
            else:
                matching_hints["nth_linear_constant_coeff_homogeneous"] = r


    # Order keys based on allhints.
    retlist = []
    for i in allhints:
        if i in matching_hints:
            retlist.append(i)


    if dict:
        # Dictionaries are ordered arbitrarily, so we need to make note of which
        # hint would come first for dsolve().  In Python 3, this should be replaced
        # with an ordered dictionary.
        matching_hints["default"] = None
        matching_hints["ordered_hints"] = tuple(retlist)
        for i in allhints:
            if i in matching_hints:
                matching_hints["default"] = i
                break
        return matching_hints
    else:
        return tuple(retlist)

@vectorize(0)
def odesimp(eq, func, order, hint):
    """
    Simplifies ODEs, including trying to solve for func and running
    constantsimp().

    It may use knowledge of the type of solution that that hint returns
    to apply additional simplifications.

    It also attempts to integrate any Integrals in the expression, if
    the hint is not an "_Integral" hint.

    This function should have no effect on expressions returned by
    dsolve(), as dsolve already calls odesimp(), but the individual hint
    functions do not call odesimp (because the dsolve() wrapper does).
    Therefore, this function is designed for mainly internal use.s

    == Example ==
        >>> from sympy import *
        >>> from sympy.solvers.ode import odesimp
        >>> x , u2, C1= symbols('x u2 C1')
        >>> f = Function('f')

        >>> eq = dsolve(x*f(x).diff(x) - f(x) - x*sin(f(x)/x), f(x),
        ... hint='1st_homogeneous_coeff_subs_indep_div_dep_Integral',
        ... simplify=False)
        >>> pprint(eq)
            x
           ----
           f(x)
             /
            |
            |                  /1 \
            |        1 + u2*sin|--|
            |                  \u2/                  /f(x)\
        -   |  -------------------------- d(u2) + log|----| = 0
            |    /          /1 \\                    \ C1 /
            |  - |1 + u2*sin|--||*u2 + u2
            |    \          \u2//
            |
           /

        >>> pprint(odesimp(eq, f(x), 1,
        ... hint='1st_homogeneous_coeff_subs_indep_div_dep'
        ... )) # (this is slow, so we skip) # doctest: +SKIP
            x
        --------- = C1
           /f(x)\
        tan|----|
           \2*x /

    """
    x = func.args[0]
    f = func.func
    C1 = Symbol('C1')

    # First, integrate, if the hint allows it.
    eq = _handle_Integral(eq, func, order, hint)

    # Second, clean up the arbitrary constants.
    # Right now, nth linear hints can put as many as 2*order constants in an
    # expression.  If that number grows with another hint, the third argument
    # here should be raise accordingly, or constantsimp() rewritten to handle
    # an arbitrary number of constants.
    eq = constantsimp(eq, x, 2*order)

    # Lastly, now that we have cleaned up the expression, try solving for func.
    # When RootOf is implemented in solve(), we will want to return a RootOf
    # everytime instead of an Equality.

    if hint[:21] == "1st_homogeneous_coeff":
        # Solutions from this hint can almost always be logcombined
        eq = logcombine(eq, assume_pos_real=True)
        if eq.lhs.is_Function and eq.lhs.func == log and eq.rhs == 0:
            eq = Eq(eq.lhs.args[0]/C1,C1)

    if eq.lhs == func and not eq.rhs.has(func):
        # The solution is already solved
        pass
    elif eq.rhs == func and not eq.lhs.has(func):
        # The solution is solved, but in reverse, so switch it
        eq = Eq(eq.rhs, eq.lhs)
    else:
        # The solution is not solved, so try to solve it
        try:
            eqsol = solve(eq, func)
            if eqsol == []:
                raise NotImplementedError
        except NotImplementedError:
            pass
        else:
            eq = [Eq(f(x), t) for t in eqsol]
            # Special handling for certain hints that we know will usually take a
            # certain form
            if hint[:21] == "1st_homogeneous_coeff":
                neweq = []
                for i in eq:
                    newi = logcombine(i, assume_pos_real=True)
                    if newi.lhs.is_Function and newi.lhs.func == log and newi.rhs == 0:
                        newi = Eq(newi.lhs.args[0]*C1,C1)
                    neweq.append(newi)
                eq = neweq
            if len(eq) == 1:
                eq = eq[0] # We only want a list if there are multiple solutions

    if hint[:25] == "nth_linear_constant_coeff":
        # Collect termms to make the solution look nice.
        # This is also necessary for constantsimp to remove unnecessary terms
        # from the particular solution from variation of parameters
        global collectterms
        sol = eq.rhs
        sol = expand_mul(sol)
        for i, reroot, imroot in collectterms:
            sol = collect(sol, x**i*exp(reroot*x)*sin(abs(imroot)*x))
            sol = collect(sol, x**i*exp(reroot*x)*cos(imroot*x))
        for i, reroot, imroot in collectterms:
            sol = collect(sol, x**i*exp(reroot*x))
        del collectterms
        eq = Eq(f(x), sol)

    # We cleaned up the costants before solving to help the solve engine with
    # a simpler expression, but the solved expression could have introduced
    # things like -C1, so rerun constantsimp() one last time before returning.
    eq = constantsimp(eq, x, 2*order)

    return eq

@vectorize(2)
def checkodesol(ode, func, sol, order='auto', solve_for_func=True):
    """
    Substitutes sol for func in ode and checks that the result is 0.

    This only works when func is one function, like f(x).  sol can be a
    single solution or a list of solutions.  Either way, each solution
    must be an Equality instance (e.g., Eq(f(x), C1*cos(x) +
    C2*sin(x))).  If it is a list of solutions, it will return a list of
    the checkodesol() result for each solution.

    It tries the following methods, in order, until it finds zero
    equivalence:

        1. Substitute the solution for f in the original equation.  This
           only works if the ode is solved for f.  It will attempt to solve
           it first unless solve_for_func == False
        2. Take n derivatives of the solution, where n is the order of
           ode, and check to see if that is equal to the solution.  This
           only works on exact odes.
        3. Take the 1st, 2nd, ..., nth derivatives of the solution, each
           time solving for the derivative of f of that order (this will
           always be possible because f is a linear operator).  Then back
           substitute each derivative into ode in reverse order.

    This function returns a tuple.  The first item in the tuple is True
    if the substitution results in 0, and False otherwise. The second
    item in the tuple is what the substitution results in.  It should
    always be 0 if the first item is True. Note that sometimes this
    function will False, but with an expression that is identically
    equal to 0, instead of returning True.  This is because simplify()
    cannot reduce the expression to 0.  If an expression returned by
    this function vanishes identically, then sol really is a solution to
    ode.

    If this function seems to hang, it is probably because of a hard
    simplification.

    To use this function to test, test the first item of the tuple.
    == Examples ==
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> f = Function('f')
        >>> checkodesol(f(x).diff(x), f(x), Eq(f(x), C1))
        (True, 0)
        >>> assert checkodesol(f(x).diff(x), f(x), Eq(f(x), C1))[0]
        >>> assert not checkodesol(f(x).diff(x), f(x), Eq(f(x), x))[0]
        >>> checkodesol(f(x).diff(x, 2), f(x), Eq(f(x), x**2))
        (False, 2)

    """
    if not isinstance(func, Function) or len(func.args) != 1:
        raise ValueError("func must be a function of one variable, not " + str(func))
    x = func.args[0]
    s = True
    testnum = 0
    if not isinstance(ode, Equality):
        ode = Eq(ode, 0)
    if not isinstance(sol, Equality):
        raise ValueError("sol must be an Equality, got " + str(sol))
    if order == 'auto':
        order = ode_order(ode, func)
    if solve_for_func and not (sol.lhs == func and not sol.rhs.has(func)) and not \
        (sol.rhs == func and not sol.lhs.has(func)):
            try:
                solved = solve(sol, func)
                if solved == []:
                    raise NotImplementedError
            except NotImplementedError:
                pass
            else:
                if len(solved) == 1:
                    result = checkodesol(ode, func, Eq(func, solved[0]), \
                        order=order, solve_for_func=False)
                else:
                    result = checkodesol(ode, func, map(lambda t: Eq(func, t), \
                        solved), order=order, solve_for_func=False)

                return result
    while s:
        if testnum == 0:
            # First pass, try substituting a solved solution directly into the ode
            # This has the highest chance of succeeding.
            if sol.lhs == func:
                    s = ode.subs(func, sol.rhs)
            elif sol.rhs == func:
                    s = ode.subs(func, sol.lhs)
            else:
                testnum += 1
                continue
            s = simplify(s.lhs - s.rhs)
            testnum += 1
        elif testnum == 1:
            # If we cannot substitute f, try seeing if the nth derivative is equal
            # This will only work for odes that are exact, by definition.
            s = simplify(trigsimp(diff(sol.lhs, x, order) - diff(sol.rhs, x, order)) - \
                trigsimp(ode.lhs) + trigsimp(ode.rhs))
 #           s2 = simplify(diff(sol.lhs, x, order) - diff(sol.rhs, x, order) - \
#                ode.lhs + ode.rhs)
            testnum += 1
        elif testnum == 2:
            # Try solving for df/dx and substituting that into the ode.
            # Thanks to Chris Smith for suggesting this method.  Many of the
            # comments below are his too.
            # The method:
            # - Take each of 1..n derivatives of the solution.
            # - Solve each nth derivative for d^(n)f/dx^(n)
            #   (the differential of that order)
            # - Back substitute into the ode in decreasing order
            #   (i.e., n, n-1, ...)
            # - Check the result for zero equivalence
            sol = sol.lhs - sol.rhs
            diffsols = {0: sol}
            for i in range(1, order + 1):
                # This is what the solution says df/dx should be.
                ds = diffsols[i - 1].diff(x)

                # Differentiation is a linear operator, so there should always
                # be 1 solution. Nonetheless, we test just to make sure.
                # We only need to solve once.  After that, we will automatically
                # have the solution to the differential in the order we want.
                if i == 1:
                    try:
                        sdf = solve(ds,func.diff(x, i))
                        if len(sdf) != 1:
                            raise NotImplementedError
                    except NotImplementedError:
                        testnum += 1
                        break
                    else:
                        diffsols[i] = sdf[0]
                else:
                    diffsols[i] = ds
            # Make sure the above didn't fail.
            if testnum > 2:
                continue
            else:
                # Substitute it into ode to check for self consistency
                for i in range(order, 0, -1):
                    # It may help if we simplify as we go
                    ode = simplify(ode.subs(func.diff(x, i), diffsols[i]))
                # No sense in overworking simplify--just prove the numerator goes to zero
                s = simplify(trigsimp((ode.lhs-ode.rhs).as_numer_denom()[0]))
                testnum += 1
        else:
            break

    if not s:
        return (True, s)
    elif s is True: # The code above never was able to change s
        raise NotImplementedError("Unable to test if " + str(sol) + \
            " is a solution to " + str(ode) + ".")
    else:
        return (False, s)

# FIXME: rewrite this as a key function, so it works in Python 3
# Python 3 removes the cmp key from sorted.  key should be better, because you
# can use it with min(), but I have no idea how to convert this into a key function
def compare_ode_sol(sol1, sol2, func, *args):
    """
    Return -1 if eq1 is simpler than eq2, 0 if they are equally complex,
    1 otherwise.

    This works like a standard Python cmp function, for use with
    functions like sort().  For example, to get the simplest expression
    from a list, you can use:

    sorted(listofodes, cmp=lambda x, y: compare_ode_sol(x, y, func))[0]

    This takes into consideration if the equations are solvable in func,
    if they contain any Integral classes (unevaluated integrals), and
    barring that, the length of the string representation of the
    expression.  Improvements to this heuristic are welcome!

    == Examples ==
        >>> from sympy import *
        >>> from sympy.solvers.ode import compare_ode_sol
        >>> x, C1 = symbols('x C1')
        >>> f = Function('f')
        >>> # # This is from dsolve(x*f(x).diff(x) - f(x) - x*sin(f(x)/x), \
        >>> # f(x), hint='1st_homogeneous_coeff_subs_indep_div_dep')
        >>> eq1 = Eq(x/tan(f(x)/(2*x)), C1)
        >>> # This is from the same ode with the
        >>> # '1st_homogeneous_coeff_subs_dep_div_indep' hint.
        >>> eq2 = Eq(x*sqrt(1 + cos(f(x)/x))/sqrt(-1 + cos(f(x)/x)), C1)
        >>> compare_ode_sol(eq1, eq2, f(x))
        -1

    """
    from sympy.core.basic import C, S

    # First, if they are the same, don't bother testing which one to use
    if sol1 == sol2:
        return 0

    # If the solutions are lists (like [Eq(f(x), sqrt(x)), Eq(f(x), -sqrt(x))],
    # then base the comparison off the worst solution in the list.
    # But when, we look at the length of the expressions at the end, use the
    # whole list.
    if isinstance(sol1, list) or isinstance(sol1, tuple):
        sol1len = sum([len(str(i)) for i in sol1])
        sol1 = sorted(sol1, cmp=lambda x, y: compare_ode_sol(x, y,
            func, *args))[len(sol1) - 1]
    else:
        sol1len = len(str(sol1))
    if isinstance(sol2, list) or isinstance(sol2, tuple):
        sol2len = sum([len(str(i)) for i in sol2])
        sol2 = sorted(sol2, cmp=lambda x, y: compare_ode_sol(x, y,
            func, *args))[len(sol2) - 1]
    else:
        sol2len = len(str(sol2))
    # Second, prefer expressions without unevaluated integrals (Integrals):
    intcmp = int(sol1.has(C.Integral)) - int(sol2.has(C.Integral))
    if intcmp:
        return intcmp

    # Next, try to solve for func.  This code will change slightly when RootOf
    # is implemented in solve().
    sol1s = 0
    sol2s = 0
    # First, see if they are already solved
    if sol1.lhs == func and not sol1.rhs.has(func) or\
        sol1.rhs == func and not sol1.lhs.has(func):
            sol1s = 1
    if sol2.lhs == func and not sol2.rhs.has(func) or\
        sol2.rhs == func and not sol2.lhs.has(func):
            sol2s = 1
    if sol2s - sol1s:
        return sol2s - sol1s
    # We are not so lucky, try solving manually
    try:
        sol1sol = solve(sol1, func)
        if sol1sol == []:
            raise NotImplementedError
    except NotImplementedError:
        pass
    else:
        sol1s = 1
    try:
        sol2sol = solve(sol2, func)
        if sol2sol == []:
            raise NotImplementedError
    except NotImplementedError:
        pass
    else:
        sol2s = 1
    if sol2s - sol1s:
        return sol2s - sol1s

    # Finally, try to return the shortest expression, naively computed
    # based on the length of the string version of the expression.  This
    # may favor combined fractions because they will not have duplicate
    # denominators, and may slightly favor expressions with fewer
    # additions and subtractions, as those are separated by spaces by
    # the printer.

    # Additional ideas for simplicity comparison are welcome, like maybe
    # checking if a equation has a larger domain.
    return cmp(sol1len, sol2len)


@vectorize(0)
def constantsimp(expr, independentsymbol, endnumber, startnumber=1,
    symbolname='C'):
    """
    Simplifies an expression with arbitrary constants in it.

    This function is written specifically to work with dsolve(), and is
    not intended for general use.

    Simplification is done by "absorbing" the arbitrary constants in to
    other arbitrary constants, numbers, and symbols that they are not
    independent of.

    The symbols must all have the same name with numbers after it, for
    example, C1, C2, C3.  The symbolname here would be 'C', the
    startnumber would be 1, and the end number would be 3.  If the
    arbitrary constants are independent of the variable x, then the
    independentsymbol would be x.  There is no need to specify the
    dependent function, such as f(x), because it already has the
    independent symbol, x, in it.

    Because terms are "absorbed" into arbitrary constants and because
    constants are renumbered after simplifying, the arbitrary constants
    in expr are not necessarily equal to the ones of the same name in
    the returned result.

    If two or more arbitrary constants are added, multiplied, or raised
    to the power of each other, they are first absorbed together into a
    single arbitrary constant.  Then the new constant is combined into
    other terms if necessary.

    Absorption is done naively.  constantsimp() does not attempt to
    expand or simplify the expression first to obtain better absorption.
    So for example, exp(C1)*exp(x) will be simplified to C1*exp(x), but
    exp(C1 + x) will be left alone.

    Constants are renumbered after simplification so that they are
    sequential, such as C1, C2, C3, and so on.  They are renumbered in
    the order that they are printed, using Basic._compare_pretty(), so
    they should be numbered in the order that they appear in an
    expression.

    In rare cases, a single constant can be "simplified" into two
    constants.  Every differential equation solution should have as many
    arbitrary constants as the order of the differential equation.  The
    result here will be technically correct, but it may, for example,
    have C1 and C2 in an expression, when C1 is actually equal to C2.
    Use your discretion in such situations, and also take advantage of
    the ability to use hints in dsolve().

    == Examples ==
        >>> from sympy import *
        >>> from sympy.solvers.ode import constantsimp
        >>> C1, C2, C3, x, y = symbols('C1 C2 C3 x y')
        >>> constantsimp(2*C1*x, x, 3)
        C1*x
        >>> constantsimp(C1 + 2 + x + y, x, 3)
        C1 + x
        >>> constantsimp(C1*C2 + 2 + x + y + C3*x, x, 3)
        C1 + x + C2*x

    """
    # We need to have an internal recursive function so that newstartnumber
    # maintains its values throughout recursive calls

    global newstartnumber
    newstartnumber = 1

    def _constantsimp(expr, independentsymbol, endnumber, startnumber=1,
    symbolname='C'):
        """
        This function works recursively.  The idea is that, for Mul,
        Add, Pow, and Function, if the class has a constant in it, then
        we can simplify it, which we do by recursing down and
        simplifying up.  Otherwise, we can skip that part of the
        expression.

        """
        from sympy.core import S
        from sympy.utilities import any
        constantsymbols = [Symbol(symbolname+"%d" % t) for t in range(startnumber,
        endnumber + 1)]
        x = independentsymbol

        if isinstance(expr, Equality):
            # For now, only treat the special case where one side of the equation
            # is a constant
            if expr.lhs in constantsymbols:
                return Eq(expr.lhs, _constantsimp(expr.rhs + expr.lhs, x, endnumber,
                startnumber, symbolname) - expr.lhs)
                # this could break if expr.lhs is absorbed into another constant,
                # but for now, the only solutions that return Eq's with a constant
                # on one side are first order.  At any rate, it will still be
                # technically correct.  The expression will just have too many
                # constants in it
            elif expr.rhs in constantsymbols:
                return Eq(_constantsimp(expr.lhs + expr.rhs, x, endnumber,
                startnumber, symbolname) - expr.rhs, expr.rhs)
            else:
                return Eq(_constantsimp(expr.lhs, x, endnumber, startnumber,
                    symbolname), _constantsimp(expr.rhs, x, endnumber,
                    startnumber, symbolname))

        if type(expr) not in (Mul, Add, Pow) and not expr.is_Function:
            # We don't know how to handle other classes
            # This also serves as the base case for the recursion
            return expr
        elif not any(expr.has(t) for t in constantsymbols):
            return expr
        else:
            newargs = []
            hasconst = False
            isPowExp = False
            reeval = False
            for i in expr.args:
                if i not in constantsymbols:
                    newargs.append(i)
                else:
                    newconst = i
                    hasconst = True
                    if expr.is_Pow and i == expr.exp:
                        isPowExp = True

            for i in range(len(newargs)):
                isimp = _constantsimp(newargs[i], x, endnumber, startnumber,
                symbolname)
                if isimp in constantsymbols:
                    reeval = True
                    hasconst = True
                    newconst = isimp
                    if expr.is_Pow and i == 1:
                        isPowExp = True
                newargs[i] = isimp
            if hasconst:
                newargs = [i for i in newargs if i.has(x)]
                if isPowExp:
                    newargs = newargs + [newconst] # Order matters in this case
                else:
                    newargs = [newconst] + newargs
            if expr.is_Pow and len(newargs) == 1:
                newargs.append(S.One)
            if expr.is_Function:
                if (len(newargs) == 0 or hasconst and len(newargs) == 1):
                    return newconst
                else:
                    newfuncargs = [_constantsimp(t, x, endnumber, startnumber,
                    symbolname) for t in expr.args]
                    return expr.new(*newfuncargs)
            else:
                newexpr = expr.new(*newargs)
                if reeval:
                    return _constantsimp(newexpr, x, endnumber, startnumber,
                    symbolname)
                else:
                    return newexpr

    def _renumber(expr, symbolname, startnumber, endnumber):
        """
        Renumber arbitrary constants in expr.

        This is a simple function that goes through and renumbers any
        Symbol with a name in the form symbolname + num where num is in
        the range from startnumber to endnumber.

        Symbols are renumbered based on Basic._compare_pretty, so they
        should be numbered roughly in the order that they appear in the
        final, printed expression.

        The structure of the function is very similar to that of
        _constantsimp().

        == Example ==
            >>> from sympy import *
            >>> from sympy.solvers.ode import constantsimp
            >>> x, C1, C2 = symbols('x C1 C2')
            >>> pprint(C2*exp(x) + C1*exp(-x))
                x       -x
            C2*e  + C1*e
            >>> pprint(constantsimp(C2*exp(x) + C1*exp(-x), x, 2))
                x       -x
            C1*e  + C2*e

        """
        from sympy.utilities import any
        constantsymbols = [Symbol(symbolname+"%d" % t) for t in range(startnumber,
        endnumber + 1)]
        global newstartnumber

        if isinstance(expr, Equality):
            return Eq(_renumber(expr.lhs, symbolname, startnumber, endnumber),
            _renumber(expr.rhs, symbolname, startnumber, endnumber))

        if type(expr) not in (Mul, Add, Pow) and not expr.is_Function and\
        not any(expr.has(t) for t in constantsymbols):
            # Base case, as above.  We better hope there aren't constants inside
            # of some other class, because they won't be simplified.
            return expr
        elif expr in constantsymbols:
            # Renumbering happens here
            newconst = Symbol(symbolname + str(newstartnumber))
            newstartnumber += 1
            return newconst
        else:
            sortedargs = list(expr.args)
            sortedargs.sort(Basic._compare_pretty)
            if expr.is_Function or expr.is_Pow:
                return expr.new(*map(lambda x: _renumber(x, symbolname, \
                startnumber, endnumber), expr.args))
            else:
                return expr.new(*map(lambda x: _renumber(x, symbolname, \
                startnumber, endnumber), sortedargs))


    simpexpr = _constantsimp(expr, independentsymbol, endnumber, startnumber,
    symbolname)

    return _renumber(simpexpr, symbolname, startnumber, endnumber)

def _handle_Integral(expr, func, order, hint):
    """
    Converts a solution with Integrals in it into an actual solution.

    For most hints, this simply runs expr.doit()

    """
    x = func.args[0]
    f = func.func
    if hint == "1st_exact":
        global exactvars
        x0 = exactvars['x0']
        y0 = exactvars['y0']
        y = exactvars['y']
        tmpsol = expr.lhs.doit()
        sol = 0
        assert tmpsol.is_Add
        for i in tmpsol.args:
            if x0 not in i and y0 not in i:
                sol += i
        assert sol != 0
        sol = Eq(sol.subs(y, f(x)),expr.rhs) # expr.rhs == C1
        del exactvars
    elif hint == "1st_exact_Integral":
        # FIXME: We still need to back substitute y
        # y = exactvars['y']
        # sol = expr.subs(y, f(x))
        # For now, we are going to have to return an expression with f(x) replaced
        # with y.  Substituting results in the y's in the second integral
        # becoming f(x), which prevents the integral from being evaluatable.
        # For example, Integral(cos(f(x)), (x, x0, x)).  If there were a way to
        # do inert substitution, that could maybe be used here instead.
        del exactvars
        sol = expr
    elif hint == "nth_linear_constant_coeff_homogeneous":
        sol = expr
    elif hint[-9:] != "_Integral":
        sol = expr.doit()
    else:
        sol = expr
    return sol

def ode_order(expr, func):
    """
    Returns the order of a given ODE with respect to func.

    This function is implemented recursively.

    == Examples ==
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> f, g = map(Function, ['f', 'g'])
        >>> ode_order(f(x).diff(x, 2) + f(x).diff(x)**2 +
        ... f(x).diff(x), f(x))
        2
        >>> ode_order(f(x).diff(x, 2) + g(x).diff(x, 3), f(x))
        2
        >>> ode_order(f(x).diff(x, 2) + g(x).diff(x, 3), g(x))
        3

    """
    a = Wild('a', exclude=[func])

    order = 0
    if isinstance(expr, Derivative) and expr.args[0] == func:
        order = len(expr.symbols)
    else:
        for arg in expr.args:
            if isinstance(arg, Derivative) and arg.args[0] == func:
                order = max(order, len(arg.symbols))
            elif expr.match(a):
                order = 0
            else :
                for arg1 in arg.args:
                    order = max(order, ode_order(arg1, func))

    return order

# FIXME: replace the general solution in the docstring with
# dsolve(equation, hint='1st_exact_Integral').  You will need to be able
# to have assumptions on P and Q that dP/dy = dQ/dx.
def ode_1st_exact(eq, func, order, match):
    r"""
    Solves 1st order exact ordinary differential equations.

    A 1st order differential equation is called exact if it is the total
    differential of a function. That is, the differential equation
    P(x, y)dx + Q(x, y)dy = 0 is exact if there is some function F(x, y)
    such that P(x, y) = dF/dx and Q(x, y) = dF/dy (d here refers to the
    partial derivative).  It can be shown that a necessary and
    sufficient condition for a first order ODE to be exact is that
    dP/dy = dQ/dx.  Then, the solution will be as given below::

        >>> from sympy import *
        >>> x, y, t, x0, y0, C1= symbols('x y t x0 y0 C1')
        >>> P, Q, F= map(Function, ['P', 'Q', 'F'])
        >>> pprint(Eq(Eq(F(x, y), Integral(P(t, y), (t, x0, x)) +
        ... Integral(Q(x0, t), (t, y0, y))), C1))
                    x                y
                    /                /
                   |                |
        F(x, y) =  |  P(t, y) dt +  |  Q(x0, t) dt = C1
                   |                |
                  /                /
                  x0               y0

    Where the first partials of P and Q exist and are continuous in a
    simply connected region.

    A note: SymPy currently has no way to represent inert substitution on
    an expression, so the hint '1st_exact_Integral' will return an integral
    with dy.  This is supposed to represent the function that you are
    solving for.

    == Example ==
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> f = Function('f')
        >>> dsolve(cos(f(x)) - (x*sin(f(x)) - f(x)**2)*f(x).diff(x),
        ... f(x), hint='1st_exact')
        x*cos(f(x)) + f(x)**3/3 == C1

    == References ==
        - http://en.wikipedia.org/wiki/Exact_differential_equation
        - M. Tenenbaum & H. Pollard, "Ordinary Differential Equations",
          Dover 1963, pp. 73

    """
    x = func.args[0]
    f = func.func
    r = match # d+e*diff(f(x),x)
    C1 = Symbol('C1')
    x0 = Symbol('x0', dummy=True)
    y0 = Symbol('y0', dummy=True)
    y = Symbol('y', dummy=True)
    global exactvars # This is the only way to pass these dummy variables to
    # _handle_Integral
    exactvars = {'y0':y0, 'x0':x0, 'y':r['y']}
    # If we ever get a Constant class, x0 and y0 should be constants, I think
    sol = C.Integral(r[r['e']].subs(x,x0),(r['y'],y0,f(x)))+C.Integral(r[r['d']],(x,x0,x))
    return Eq(sol, C1)


def ode_1st_homogeneous_coeff_best(eq, func, order, match):
    r"""
    Returns the best solution to an ODE from the two hints
    '1st_homogeneous_coeff_subs_dep_div_indep' and
    '1st_homogeneous_coeff_subs_indep_div_dep'.

    This is as determined by compare_ode_sol().

    See the ode_1st_homogeneous_coeff_subs_indep_div_dep() and
    ode_1st_homogeneous_coeff_subs_dep_div_indep() docstrings for more
    information on these hints.  Note that there is no
    '1st_homogeneous_coeff_best_Integral' hint.

    == Example ==
    ::
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> f = Function('f')
        >>> pprint(dsolve(2*x*f(x) + (x**2 + f(x)**2)*f(x).diff(x), f(x),
        ... hint='1st_homogeneous_coeff_best'))
              ___________
             /         2
            /       3*x
           /   1 + ----- *f(x) = C1
        3 /         2
        \/         f (x)

    == References ==
        - http://en.wikipedia.org/wiki/Homogeneous_differential_equation
        - M. Tenenbaum & H. Pollard, "Ordinary Differential Equations",
          Dover 1963, pp. 59

    """
    # There are two substitutions that solve the equation, u1=y/x and u2=x/y
    # They produce different integrals, so try them both and see which
    # one is easier.
    sol1 = ode_1st_homogeneous_coeff_subs_indep_div_dep(eq,
    func, order, match)
    sol2 = ode_1st_homogeneous_coeff_subs_dep_div_indep(eq,
    func, order, match)
    if match.get('simplify', True):
        sol1 = odesimp(sol1, func, order, "1st_homogeneous_coeff_subs_indep_div_dep")
        sol2 = odesimp(sol2, func, order, "1st_homogeneous_coeff_subs_dep_div_indep")
    return sorted([sol1, sol2], cmp=lambda x, y: compare_ode_sol(x, y, func))[0]

def ode_1st_homogeneous_coeff_subs_dep_div_indep(eq, func, order, match):
    r"""
    Solves a 1st order differential equation with homogeneous coefficients
    using the substitution
    u1 = <dependent variable>/<independent variable>.

    This is a differential equation P(x, y) + Q(x, y)dy/dx = 0, that P
    and Q are homogeneous of the same order.  A function F(x, y) is
    homogeneous of order n if F(xt, yt) = t**n*F(x, y).  Equivalently,
    F(x, y) can be rewritten as G(y/x) or H(x/y).  See also the
    docstring of homogeneous_order().

    If the coefficients P and Q in the  differential equation above are
    homogeneous functions of the same order, then it can be shown that
    the substitution y = u1*x (u1 = y/x) will turn the differential
    equation into an equation separable in the variables x and u.  if
    h(u1) is the function that results from making the substitution
    u1 = f(x)/x on P(x, f(x)) and g(u2) is the function that results
    from the substitution on Q(x, f(x)) in the differential equation
    P(x, f(x)) + Q(x, f(x))*diff(f(x), x) = 0, then the general solution
    is::

        >>> from sympy import *
        >>> x = Symbol('x')
        >>> f, g, h = map(Function, ['f', 'g', 'h'])
        >>> pprint(dsolve(g(f(x)/x) + h(f(x)/x)*f(x).diff(x), f(x),
        ... hint='1st_homogeneous_coeff_subs_dep_div_indep_Integral'))
           f(x)
           ----
            x
             /
            |
            |       -h(u1)
        -   |  ---------------- d(u1) + log(C1*x) = 0
            |  u1*h(u1) + g(u1)
            |
           /


    Where u1*h(u1) + g(u1) != 0 and x != 0.

    See also the docstrings of ode_1st_homogeneous_coeff_best() and
    ode_1st_homogeneous_coeff_subs_indep_div_dep().

    == Example ==
    ::
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> f = Function('f')
        >>> pprint(dsolve(2*x*f(x) + (x**2 + f(x)**2)*f(x).diff(x), f(x),
        ... hint='1st_homogeneous_coeff_subs_dep_div_indep'))
                ________________
               /           3
              /  3*f(x)   f (x)
        x*   /   ------ + -----  = C1
          3 /      x         3
          \/                x

    == References ==
        - http://en.wikipedia.org/wiki/Homogeneous_differential_equation
        - M. Tenenbaum & H. Pollard, "Ordinary Differential Equations",
          Dover 1963, pp. 59

    """
    x = func.args[0]
    f = func.func
    u1 = Symbol('u1', dummy=True) # u1 == f(x)/x
    r = match # d+e*diff(f(x),x)
    C1 = Symbol('C1')
    int = C.Integral((-r[r['e']]/(r[r['d']]+u1*r[r['e']])).subs({x:1, r['y']:u1}),
        (u1, None, f(x)/x))
    sol = logcombine(Eq(log(x), int + log(C1)), assume_pos_real=True)
    return sol

def ode_1st_homogeneous_coeff_subs_indep_div_dep(eq, func, order, match):
    r"""
    Solves a 1st order differential equation with homogeneous coefficients
    using the substitution
    u2 = <independent variable>/<dependent variable>.

    This is a differential equation P(x, y) + Q(x, y)dy/dx = 0, that P
    and Q are homogeneous of the same order.  A function F(x, y) is
    homogeneous of order n if F(xt, yt) = t**n*F(x, y).  Equivalently,
    F(x, y) can be rewritten as G(y/x) or H(x/y).  See also the
    docstring of homogeneous_order().

    If the coefficients P and Q in the  differential equation above are
    homogeneous functions of the same order, then it can be shown that
    the substitution x = u2*y (u2 = x/y) will turn the differential
    equation into an equation separable in the variables y and u2.  if
    h(u2) is the function that results from making the substitution
    u2 = x/f(x) on P(x, f(x)) and g(u2) is the function that results
    from the substitution on Q(x, f(x)) in the differential equation
    P(x, f(x)) + Q(x, f(x))*diff(f(x), x) = 0, then the general solution
    is:

    >>> from sympy import *
    >>> x = Symbol('x')
    >>> f, g, h = map(Function, ['f', 'g', 'h'])
    >>> pprint(dsolve(g(x/f(x)) + h(x/f(x))*f(x).diff(x), f(x),
    ... hint='1st_homogeneous_coeff_subs_indep_div_dep_Integral'))
                 x
                ----
                f(x)
                  /
                 |
                 |       -g(u2)
                 |  ---------------- d(u2)
                 |  u2*g(u2) + h(u2)
                 |
                /
    f(x) = C1*e


    Where u2*g(u2) + h(u2) != 0 and f(x) != 0.

    See also the docstrings of ode_1st_homogeneous_coeff_best() and
    ode_1st_homogeneous_coeff_subs_dep_div_indep().

    == Example ==
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> f = Function('f')
        >>> pprint(dsolve(2*x*f(x) + (x**2 + f(x)**2)*f(x).diff(x), f(x),
        ... hint='1st_homogeneous_coeff_subs_indep_div_dep'))
              ___________
             /         2
            /       3*x
           /   1 + ----- *f(x) = C1
        3 /         2
        \/         f (x)

    == References ==
        - http://en.wikipedia.org/wiki/Homogeneous_differential_equation
        - M. Tenenbaum & H. Pollard, "Ordinary Differential Equations",
          Dover 1963, pp. 59

    """
    x = func.args[0]
    f = func.func
    u2 = Symbol('u2', dummy=True) # u2 == x/f(x)
    r = match # d+e*diff(f(x),x)
    C1 = Symbol('C1')
    int = C.Integral((-r[r['d']]/(r[r['e']]+u2*r[r['d']])).subs({x:u2, r['y']:1}),
        (u2, None, x/f(x)))
    sol = logcombine(Eq(log(f(x)), int + log(C1)), assume_pos_real=True)
    return sol

# XXX: Should this function maybe go somewhere else?
def homogeneous_order(eq, *symbols):
    """
    Returns the order n if g is homogeneous and None if it is not
    homogeneous.

    Determines if a function is homogeneous and if so of what order.
    A function f(x,y,...) is homogeneous of order n if
    f(t*x,t*y,t*...) == t**n*f(x,y,...).  The function is implemented recursively.

    If the function is of two variables, F(x, y), then f being
    homogeneous of any order is equivalent to being able to rewrite
    F(x, y) as G(x/y) or H(y/x).  This fact is used to solve 1st order
    ordinary differential equations whose coefficients are homogeneous
    of the same order (see the docstrings of
    ode.ode_1st_homogeneous_coeff_subs_indep_div_dep() and
    ode.ode_1st_homogeneous_coeff_subs_indep_div_dep()

    Symbols can be functions, but every argument of the function must be
    a symbol, and the arguments of the function that appear in the
    expression must match those given in the list of symbols.  If a
    declared function appears with different arguments than given in the
    list of symbols, None is returned.

    Example::
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> y = Symbol('y')
        >>> f = Function('f')
        >>> homogeneous_order(f(x), f(x)) == None
        True
        >>> homogeneous_order(f(x,y), f(y, x), x, y) == None
        True
        >>> homogeneous_order(f(x), f(x), x)
        1

    == Examples ==
        >>> homogeneous_order(x**2*f(x)/sqrt(x**2+f(x)**2), x, f(x))
        2
        >>> homogeneous_order(x**2+f(x), x, f(x)) == None
        True

    """
    if eq.has(log):
        eq = logcombine(eq, assume_pos_real=True)
    return _homogeneous_order(eq, *symbols)

def _homogeneous_order(eq, *symbols):
    """
    The real work for homogeneous_order.

    This runs as a separate function call so that logcombine doesn't
    endlessly put back together what homogeneous_order is trying to take
    apart.
    """
    from sympy.utilities import all, any
    if not symbols:
        raise ValueError, "homogeneous_order: no symbols were given."

    n = set()

    # Replace all functions with dummy variables

    if any(getattr(i, 'is_Function') for i in symbols):
        for i in symbols:
            if i.is_Function:
                if not all(map((lambda i: i in symbols), i.args)):
                    return None
                elif i not in symbols:
                    pass
                else:
                    dummyvar = numbered_symbols(prefix='d', dummy=True).next()
                    eq = eq.subs(i, dummyvar)
                    symbols = list(symbols)
                    symbols.remove(i)
                    symbols.append(dummyvar)
                    symbols = tuple(symbols)

    # The following are not supported
    if eq.is_Order or eq.is_Derivative:
        return None

    # These are all constants
    if type(eq) in (int, float) or eq.is_Number or eq.is_Integer or \
    eq.is_Rational or eq.is_NumberSymbol or eq.is_Real:
        return sympify(0)

    # Break the equation into additive parts
    if eq.is_Add:
        s = set()
        for i in eq.args:
            s.add(_homogeneous_order(i, *symbols))
        if len(s) != 1:
            return None
        else:
            n = s

    if eq.is_Pow:
        if not eq.args[1].is_Number:
            return None
        o = _homogeneous_order(eq.args[0], *symbols)
        if o == None:
            return None
        else:
            n.add(sympify(o*eq.args[1]))

    t = Symbol('t', dummy=True, positive=True) # It is sufficient that t > 0
    r = Wild('r', exclude=[t])
    a = Wild('a', exclude=[t])
    eqs = eq.subs(dict(zip(symbols,(t*i for i in symbols))))

    if eqs.is_Mul:
        if t not in eqs:
            n.add(sympify(0))
        else:
            m = eqs.match(r*t**a)
            if m:
                n.add(sympify(m[a]))
            else:
                s = 0
                for i in eq.args:
                    o = _homogeneous_order(i, *symbols)
                    if o == None:
                        return None
                    else:
                        s += o
                n.add(sympify(s))

    if eq.is_Function:
        if eq.func == log:
            # The only possibility to pull a t out of a function is a power in
            # a logarithm.  This is very likely due to calling of logcombine().
            if eq.args[0].is_Pow:
                return _homogeneous_order(eq.args[0].args[1]*log(eq.args[0].args[0]),\
                    *symbols)
            elif eq.args[0].is_Mul and all(i.is_Pow for i in iter(eq.args[0].args)):
                arg = 1
                pows = set()
                for i in eq.args[0].args:
                    if i.args[1].args[0] == -1:
                        arg *= 1/i.args[0]
                        pows.add(sympify(-1*i.args[1]))
                    else:
                        arg *= i.args[0]
                        pows.add(sympify(i.args[1]))
                if len(pows) != 1:
                    return None
                else:
                    return _homogeneous_order(pows.pop()*log(arg), *symbols)
            else:
                if _homogeneous_order(eq.args[0], *symbols) == 0:
                    return sympify(0)
                else:
                    return None
        else:
            if _homogeneous_order(eq.args[0], *symbols) == 0:
                return sympify(0)
            else:
                return None

    if len(n) != 1 or n == None:
        return None
    else:
        return n.pop()

    return None



def ode_1st_linear(eq, func, order, match):
    r"""
    Solves 1st order linear differential equations.

    These are differential equations of the form dy/dx _ P(x)*y = Q(x).
    These kinds of differential equations can be solved in a general
    way.  The integrating factor exp(Integral(P(x), x)) will turn the
    equation into a separable equation.  The general solution is::

        >>> from sympy import *
        >>> x = Symbol('x')
        >>> f, P, Q = map(Function, ['f', 'P', 'Q'])
        >>> pprint(dsolve(Eq(f(x).diff(x) + P(x)*f(x), Q(x)), f(x),
        ... hint='1st_linear_Integral'))
               /       /                   \
               |      |                    |
               |      |         /          |     /
               |      |        |           |    |
               |      |        | P(x) dx   |  - | P(x) dx
               |      |        |           |    |
               |      |       /            |   /
        f(x) = |C1 +  | Q(x)*e           dx|*e
               |      |                    |
               \     /                     /


    == Example ==
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> f = Function('f')
        >>> pprint(dsolve(Eq(x*diff(f(x), x) - f(x), x**2*sin(x)),
        ... f(x), '1st_linear'))
        f(x) = x*(C1 - cos(x))

    == References ==
        - http://en.wikipedia.org/wiki/Linear_differential_equation#First_order_equation
        - M. Tenenbaum & H. Pollard, "Ordinary Differential Equations",
          Dover 1963, pp. 92

    """
    x = func.args[0]
    f = func.func
    r = match # a*diff(f(x),x) + b*f(x) + c
    C1 = Symbol('C1')
    t = exp(C.Integral(r[r['b']]/r[r['a']], x))
    tt = C.Integral(t*(-r[r['c']]/r[r['a']]), x)
    return Eq(f(x),(tt + C1)/t)

def ode_Bernoulli(eq, func, order, match):
    r"""
    Solves Bernoulli differential equations.

    These are equations of the form dy/dx + P(x)*y = Q(x)*y**n, n != 1.
    The substitution w = 1/y**(1-n) will transform an equation of this
    form into one that is linear (see the docstring of
    ode_1st_linear()).  The general solution is::

        >>> from sympy import *
        >>> x, n = symbols('x n')
        >>> f, P, Q = map(Function, ['f', 'P', 'Q'])
        >>> pprint(dsolve(Eq(f(x).diff(x) + P(x)*f(x), Q(x)*f(x)**n),
        ... f(x), hint='Bernoulli_Integral')) # doctest: +SKIP
                                                                                       1
                                                                                      ----
                                                                                     1 - n
               //                /                            \                     \
               ||               |                             |                     |
               ||               |                  /          |             /       |
               ||               |                 |           |            |        |
               ||               |        (1 - n)* | P(x) dx   |  (-1 + n)* | P(x) dx|
               ||               |                 |           |            |        |
               ||               |                /            |           /         |
        f(x) = ||C1 + (-1 + n)* | -Q(x)*e                   dx|*e                   |
               ||               |                             |                     |
               \              /                              /                     /


    Note that when n = 1, then the equation is separable (see the
    docstring of ode_separable()).

    >>> pprint(dsolve(Eq(f(x).diff(x) + P(x)*f(x), Q(x)*f(x)), f(x),
    ... hint='separable_Integral'))
     f(x)
       /
      |                /
      |  1            |
      |  - dy = C1 +  | (-P(x) + Q(x)) dx
      |  y            |
      |              /
     /


    == Example ==
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> f = Function('f')
        >>> pprint(dsolve(Eq(x*f(x).diff(x) + f(x), log(x)*f(x)**2),
        ... f(x), hint='Bernoulli')) # doctest: +SKIP
                        1
        f(x) = -------------------
                 /     log(x)   1\
               x*|C1 + ------ + -|
                 \       x      x/

    == References ==
        - http://en.wikipedia.org/wiki/Bernoulli_differential_equation
        - M. Tenenbaum & H. Pollard, "Ordinary Differential Equations",
          Dover 1963, pp. 95

    """
    x = func.args[0]
    f = func.func
    r = match # a*diff(f(x),x) + b*f(x) + c*f(x)**n, n != 1
    C1 = Symbol('C1')
    t = exp((1-r[r['n']])*C.Integral(r[r['b']]/r[r['a']],x))
    tt = (r[r['n']]-1)*C.Integral(t*r[r['c']]/r[r['a']],x)
    return Eq(f(x),((tt + C1)/t)**(1/(1-r[r['n']])))

def ode_Liouville(eq, func, order, match):
    r"""
    Solves 2nd order Liouville differential equations.

    The general form of a Liouville ODE is
    d^2y/dx^2 + g(y)*(dy/dx)**2 + h(x)*dy/dx.  The general solution is::
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> f, g, h = map(Function, ['f', 'g', 'h'])
        >>> pprint(dsolve(Eq(diff(f(x),x,x) + g(f(x))*diff(f(x),x)**2 +
        ... h(x)*diff(f(x),x), 0), f(x), hint='Liouville_Integral'))
                                          f(x)
                  /                     /
                 |                     |
                 |     /               |     /
                 |    |                |    |
                 |  - | h(x) dx        |    | g(y) dy
                 |    |                |    |
                 |   /                 |   /
        C1 + C2* | e            dx +   |  e           dy = 0
                 |                     |
                /                     /

    == Example ==
    ::
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> f = Function('f')
        >>> pprint(dsolve(diff(f(x), x, x) + diff(f(x), x)**2/f(x) +
        ... diff(f(x), x)/x, f(x), hint='Liouville'))
                   ________________           ________________
        [f(x) = -\/ C1 + C2*log(x) , f(x) = \/ C1 + C2*log(x) ]


    == References ==
        - Goldstein and Braun, "Advanced Methods for the Solution of
          Differential Equations", pp. 98
        - http://www.maplesoft.com/support/help/view.aspx?path=odeadvisor/Liouville

    """
    # Liouville ODE f(x).diff(x, 2) + g(f(x))*(f(x).diff(x, 2))**2 + h(x)*f(x).diff(x)
    # See Goldstein and Braun, "Advanced Methods for the Solution of
    # Differential Equations", pg. 98, as well as
    # http://www.maplesoft.com/support/help/view.aspx?path=odeadvisor/Liouville
    x = func.args[0]
    f = func.func
    r = match # f(x).diff(x, 2) + g*f(x).diff(x)**2 + h*f(x).diff(x)
    y = r['y']
    C1 = Symbol('C1')
    C2 = Symbol('C2')
    int = C.Integral(exp(C.Integral(r['g'], y)), (y, None, f(x)))
    sol = Eq(int + C1*C.Integral(exp(-C.Integral(r['h'], x)), x) + C2, 0)
    return sol


def _nth_linear_match(eq, func, order):
    """
    Matches a differential equation to the linear form:

    a_n(x)y^(n) + ... + a_1(x)y' + a_0(x)y + B(x) = 0

    Returns a dict of order:coeff terms, where order is the order of the
    derivative on each term, and coeff is the coefficient of that
    derivative.  The key 'b' holds the function B(x). Returns None if
    the ode is not linear.  This function assumes that func has already
    been checked to be good.

    == Examples ==
        >>> from sympy.solvers.ode import _nth_linear_match
        >>> x = Symbol('x')
        >>> f = Function('f')
        >>> _nth_linear_match(f(x).diff(x, 3) + 2*f(x).diff(x) +
        ... x*f(x).diff(x, 2) + cos(x)*f(x).diff(x) + x - f(x) -
        ... sin(x), f(x), 3)
        {'b': x - sin(x), 1: 2 + cos(x), 0: -1, 2: x, 3: 1}
        >>> _nth_linear_match(f(x).diff(x, 3) + 2*f(x).diff(x) +
        ... x*f(x).diff(x, 2) + cos(x)*f(x).diff(x) + x - f(x) -
        ... sin(f(x)), f(x), 3) == None
        True

    """
    from sympy import S
    x = func.args[0]
    terms={'b': S.Zero}
    for i in range(order + 1):
        terms[i] = S.Zero
    # FIXME: use .as_Add() here, when the new polys module is merged in.
    if eq.is_Add:
        eqargs = eq.args
    else:
        eqargs = [eq]
    for i in eqargs:
        if not i.has(func):
            terms['b'] += i
        else:
             # .coeff(func) gets func and all derivatives of func
            c = i.coeff(func)
            if not c:
                return None
            else:
                t = i.extract_multiplicatively(c)
                if t == func:
                    terms[0] += c
                elif isinstance(t, Derivative):
                    # Make sure every symbol is x
                    if not all(map(lambda t: t == x, t.symbols)) or \
                        not t.expr == func:
                            return None
                    else:
                        terms[len(t.symbols)] += c
                else:
                    return None
    return terms

def ode_nth_linear_constant_coeff_homogeneous(eq, func, order, match, returns='sol'):
    """
    Solves an nth order linear homogeneous differential equation with
    constant coefficients.

    This is an equation of the form a_n*f(x)^(n) + a_(n-1)*f(x)^(n-1) +
    ... + a1*f'(x) + a0*f(x) = 0

    These equations can be solved in a general mannar, by taking the
    roots of the characteristic equation a_n*m**n + a_(n-1)*m**(n-1) +
    ... + a1*m + a0 = 0.  The solution will then be the sum of
    Cn*x**i*exp(r*x) terms, for each  where Cn is an arbitrary constant,
    r is a root of the characteristic equation and i is is one of each
    from 0 to the multiplicity of the root - 1 (for example, a root 3 of
    multilpicity 2 would create the terms C1*exp(3*x) + C2*x*exp(3*x)).
    The exponential is usually expanded for complex roots using Euler's
    equation exp(I*x) = cos(x) + I*sin(x).  Complex roots always come in
    conjugate pars in polynomials with real coefficients, so the two
    roots will be represented (after simplifying the constants) as
    exp(a*x)*(C1*cos(b*x) + C2*sin(b*x)).

    If SymPy cannot find exact roots to the characteristic equation, a
    RootOf instance will be return in it's stead.

    >>> from sympy import *
    >>> x = Symbol('x')
    >>> f = Function('f')
    >>> dsolve(f(x).diff(x, 5) + 10*f(x).diff(x) - 2*f(x), f(x),
    ... hint='nth_linear_constant_coeff_homogeneous')
    f(x) == C1*exp(x*RootOf(_m**5 + 10*_m - 2, _m, index=0)) +
    C2*exp(x*RootOf(_m**5 + 10*_m - 2, _m, index=1)) +
    C3*exp(x*RootOf(_m**5 + 10*_m - 2, _m, index=2)) +
    C4*exp(x*RootOf(_m**5 + 10*_m - 2, _m, index=3)) +
    C5*exp(x*RootOf(_m**5 + 10*_m - 2, _m, index=4))

    Note that because this method does not involve integration, there is
    no 'nth_linear_constant_coeff_homogeneous_Integral' hint.

    The following is for internal use:

    returns = 'sol' returns the solution to the ODE.
    returns = 'list'  returns a list of linearly independent solutions,
    for use with non homogeneous solution methods like variation of
    parameters and undetermined coefficients.  Note that, though the
    solutions should be linearly independent, this function does not
    explicitly check that.  You can do "assert
    simplify(wronskian(sollist)) != 0" to check for linear independence.
     Also, "assert len(sollist) == order" will need to pass.
     returns = 'both', return a dictionary {'sol':solution to ODE,
     'list': list of linearly independent solutions}.

    == Example ==
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> f = Function('f')
        >>> pprint(dsolve(f(x).diff(x, 4) + 2*f(x).diff(x, 3) -
        ... 2*f(x).diff(x, 2) - 6*f(x).diff(x) + 5*f(x), f(x),
        ... hint='nth_linear_constant_coeff_homogeneous'))
                            x                            -2*x
        f(x) = (C1 + C2*x)*e  + (C3*cos(x) + C4*sin(x))*e


    == References ==
        - http://en.wikipedia.org/wiki/Linear_differential_equation
            section: Nonhomogeneous_equation_with_constant_coefficients
        - M. Tenenbaum & H. Pollard, "Ordinary Differential Equations",
          Dover 1963, pp. 211

    """
    from sympy.core.basic import S
    x = func.args[0]
    f = func.func
    r = match
    # A generator of constants
    constants = numbered_symbols(prefix='C', function=Symbol, start=1)
    # First, set up characteristic equation.
    m = Symbol('m', dummy=True)
    chareq = S.Zero
    for i in r.keys():
        if type(i) == str:
            pass
        else:
            chareq += r[i]*m**i
    chareqroots = RootsOf(chareq, m)
    charroots_exact = list(chareqroots.exact_roots())
    charroots_formal = list(chareqroots.formal_roots())
    if charroots_formal and discriminant(chareq, m) == 0:
        # If Poly cannot find the roots explicitly, we can only return
        # an expression in terms of RootOf's if we know the roots
        # are not repeated.  We use the fact that a polynomial has
        # repeated roots iff its discriminant == 0.

        # Ideally, RootOf would cancel out roots from charroots_exact, so
        # we check the discriminant of only the unknown part of the chareq.
        # See issue 1557.
        raise NotImplementedError("Cannot find all of the roots of " + \
        "characteristic equation " + str(chareq) + ", which has " + \
        "repeated roots.")
    # Create a dict root: multiplicity or charroots
    charroots = {}
    for i in charroots_exact + charroots_formal:
        if i in charroots:
            charroots[i] += 1
        else:
            charroots[i] = 1
    gsol = S(0)
    # We need keep track of terms so we can run collect() at the end.
    # This is necessary for constantsimp to work properly.
    global collectterms
    collectterms = []
    for root, multiplicity in charroots.items():
        for i in range(multiplicity):
            if isinstance(root, RootOf):
                gsol += exp(root*x)*constants.next()
                assert multiplicity == 1
                collectterms = [(0, root, 0)] + collectterms
            else:
                reroot = re(root)
                imroot = im(root)
                gsol += x**i*exp(reroot*x)*(constants.next()*sin(abs(imroot)*x) \
                + constants.next()*cos(imroot*x))
                # This ordering is important
                collectterms = [(i, reroot, imroot)] + collectterms
    if returns == 'sol':
        return Eq(f(x), gsol)
    elif returns == 'list' or returns == 'both':
        # Create a list of (hopefully) linearly independent solutions
        gensols = []
        # Keep track of when to use sin or cos for nonzero imroot
        trigdict = {}
        for i, reroot, imroot in collectterms:
            if imroot == 0:
                gensols.append(x**i*exp(reroot*x))
            else:
                if x**i*exp(reroot*x)*sin(abs(imroot)*x) in gensols:
                    gensols.append(x**i*exp(reroot*x)*cos(imroot*x))
                else:
                    gensols.append(x**i*exp(reroot*x)*sin(abs(imroot)*x))
        if returns == 'list':
            return gensols
        else:
            return {'sol':Eq(f(x), gsol), 'list':gensols}
    else:
        raise ValueError('Unknown value for key "returns".')

def ode_nth_linear_constant_coeff_undetermined_coefficients(eq, func, order, match):
    r"""
    Solves an nth order linear differential equation with constant
    coefficients using the method of undetermined coefficients.

    This method works on differential equations of the form a_n*f(x)^(n)
    + a_(n-1)*f(x)^(n-1) + ... + a1*f'(x) + a0*f(x) = P(x), where P(x)
    is a function that has a finite number of linearly independent
    derivatives.

    Functions that fit this requirement are finite sums functions of the
    form a*x**i*exp(b*x)*sin(c*x + d) or a*x**i*exp(b*x)*cos(c*x + d),
    where i is a non-negative integer and a, b, c, and d are constants.
    For example any polynomial in x, functions like x**2*exp(2*x),
    x*sin(x), and exp(x)*cos(x) can all be used.  Products of sin's and
    cos's have a finite number of derivatives, because they can be
    expanded into sin(a*x) and cos(b*x) terms.  However, SymPy currently
    cannot do that expansion, so you will need to manually rewrite the
    expression in terms of the above to use this method.  So, for example,
    you will need to manually convert sin(x)**2 into (1 + cos(2*x))/2 to
    properly apply the method of undetermined coefficients on it.

    This method works by creating a trial function from the expression
    and all of its linear independent derivatives and substituting them
    into the original ODE.  The coefficients for each term will be a
    system of linear equations, which are be solved for and substituted,
    giving the solution.  If any of the trial functions are linearly
    dependent on the solution to the homogeneous equation, they are
    multiplied by sufficient x to make them linearly independent.

    == Example ==
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> f = Function('f')
        >>> pprint(dsolve(f(x).diff(x, 2) + 2*f(x).diff(x) + f(x) -
        ... 4*exp(-x)*x**2 + cos(2*x), f(x),
        ... hint='nth_linear_constant_coeff_undetermined_coefficients'))
                                           /             4\
                 4*sin(2*x)   3*cos(2*x)   |            x |  -x
        f(x) = - ---------- + ---------- + |C1 + C2*x + --|*e
                     25           25       \            3 /

    == References ==
        - http://en.wikipedia.org/wiki/Method_of_undetermined_coefficients
        - M. Tenenbaum & H. Pollard, "Ordinary Differential Equations",
          Dover 1963, pp. 221

    """
    gensol = ode_nth_linear_constant_coeff_homogeneous(eq, func, order, match,
        returns='both')
    match.update(gensol)
    return _solve_undetermined_coefficients(eq, func, order, match)

def _solve_undetermined_coefficients(eq, func, order, match):
    """
    Helper function for the method of undetermined coefficients.

    See the ode_nth_linear_constant_coeff_undetermined_coefficients()
    docstring for more information on this method.

    match should be a dictionary that has the following keys:
    'list' - A list of solutions to the homogeneous equation, such as
         the list returned by
         ode_nth_linear_constant_coeff_homogeneous(returns='list')
    'sol' - The general solution, such as the solution returned by
        ode_nth_linear_constant_coeff_homogeneous(returns='sol')
    'trialset' - The set of trial functions as returned by
        _undetermined_coefficients_match()['trialset']

    """
    x = func.args[0]
    f = func.func
    r = match
    coeffs = numbered_symbols('a', dummy=True)
    coefflist = []
    gensols = r['list']
    gsol = r['sol']
    trialset = r['trialset']
    notneedset = set([])
    newtrialset = set([])
    global collectterms
    if len(gensols) != order:
        raise NotImplementedError("Cannot find " + str(order) + \
        " solutions to the homogeneous equation nessesary to apply " + \
        "undetermined coefficients to " + str(eq) + "(number of terms != order)")
    usedsin = set([])
    mult = 0 # The multiplicity of the root
    getmult = True
    for i, reroot, imroot in collectterms:
        if getmult:
            mult = i + 1
            getmult = False
        if i == 0:
            getmult = True
        if imroot:
            # Alternate between sin and cos
            if (i, reroot) in usedsin:
                check = x**i*exp(reroot*x)*cos(imroot*x)
            else:
                check = x**i*exp(reroot*x)*sin(abs(imroot)*x)
                usedsin.add((i, reroot))
        else:
            check = x**i*exp(reroot*x)

        if check in trialset:
            # If an element of the trial function is already part of the homogeneous
            # solution, we need to multiply by sufficient x to make it linearly
            # independent.  We also don't need to bother checking for the coefficients
            # on those elements, since we already know it will be 0.
            while True:
                if check*x**mult in trialset:
                    mult += 1
                else:
                    break
            trialset.add(check*x**mult)
            notneedset.add(check)

    newtrialset = trialset - notneedset

    trialfunc = 0
    for i in newtrialset:
        c = coeffs.next()
        coefflist.append(c)
        trialfunc += c*i

    eqs = eq.subs(f(x), trialfunc)
    coeffsdict = dict(zip(trialset, [0 for i in range(len(trialset) + 1)]))
    if eqs.is_Add:
        eqs = expand_mul(eqs)
        for i in eqs.args:
            s = separatevars(i, dict=True, symbols=[x])
            if coeffsdict.has_key(s[x]):
                coeffsdict[s[x]] += s['coeff']
            else:
                # We removed that term above because we already know its coeff
                # will be 0
                pass
    else:
        s = separatevars(eqs, dict=True, symbols=[x])
        coeffsdict[s[x]] += s['coeff']

    coeffvals = solve(coeffsdict.values(), coefflist)

    if not coeffvals:
        raise NotImplementedError("Could not solve " + str(eq) + " using the " + \
            " method of undetermined coefficients (unable to solve for coefficients).")

    psol = trialfunc.subs(coeffvals)

    return Eq(f(x), gsol.rhs + psol)

def _undetermined_coefficients_match(expr, x):
    """
    Returns a trial function match if undetermined coefficients can be
    applied to expr, and None otherwise.

    A trial expression can be found for an expression for use with the
    method of undetermined coefficients if the expression is an
    additive/multiplicative combination of constants, polynomials in x
    (the independent variable of expr), sin(a*x + b), cos(a*x + b), and
    exp(a*x) terms (in other words, it has a finite number of linearly
    independent derivatives).

    Note that you may still need to multiply each term returned here by
    sufficient x to make it linearly independent with the solutions to
    the homogeneous equation.

    This is intended for internal use by undetermined_coefficients
    hints.

    SymPy currently has no way to convert sin(x)**n*cos(y)**m into a sum
    of only sin(a*x) and cos(b*x) terms, so these are not implemented.
    So, for example, you will need to manually convert sin(x)**2 into
    (1 + cos(2*x))/2 to properly apply the method of undetermined
    coefficients on it.

    == Example ==
        >>> from sympy import *
        >>> from sympy.solvers.ode import _undetermined_coefficients_match
        >>> x = Symbol('x')
        >>> _undetermined_coefficients_match(9*x*exp(x) + exp(-x), x)
        {'test': True, 'trialset': set([x*exp(x), exp(x), exp(-x)])}
        >>> _undetermined_coefficients_match(log(x), x)
        {'test': False}

    """
    from sympy import S
    a = Wild('a', exclude=[x])
    b = Wild('b', exclude=[x])
    expr = powsimp(expr, combine='exp') # exp(x)*exp(2*x + 1) => exp(3*x + 1)
    retdict = {}
    def _test_term(expr, x):
        """
        Test if expr fits the proper form for undetermined coefficients.
        """
        if expr.is_Add:
            return all([_test_term(i, x) for i in expr.args])
        elif expr.is_Mul:
            if expr.has(sin) or expr.has(cos):
                foundtrig = False
                # Make sure that there is only on trig function in the args.
                # See the docstring.
                for i in expr.args:
                    if i.has(sin) or i.has(cos):
                        if foundtrig:
                            return False
                        else:
                            foundtrig = True
            return all([_test_term(i, x) for i in expr.args])
        elif expr.is_Function:
            if expr.func in (sin, cos, exp):
                if expr.args[0].match(a*x + b):
                    return True
                else:
                    return False
            else:
                return False
        elif expr.is_Pow and expr.base.is_Symbol and expr.exp.is_Integer and \
            expr.exp >= 0:
                return True
        elif expr.is_Pow and expr.base.is_number:
            if expr.exp.match(a*x + b):
                return True
            else:
                return False
        elif expr.is_Symbol or expr.is_Number:
            return True
        else:
            return False

    def _get_trial_set(expr, x, exprs=set([])):
        """
        Returns a set of trial terms for undetermined coefficients.

        The idea behind undetermined coefficients is that the terms
        expression repeat themselves after a finite number of
        derivatives, except for the coefficients (they are linearly
        dependent).  So if we collect these, we should have the terms of
        our trial function.
        """
        def _remove_coefficient(expr, x):
            """
            Returns the expression without a coefficient.

            Similar to expr.as_independent(x)[1], except it only works
            multiplicatively.
            """
            # I was using the below match, but it doesn't always put all of the
            # coefficient in c.  c.f. 2**x*6*exp(x)*log(2)
            # The below code is probably cleaner anyway.
#            c = Wild('c', exclude=[x])
#            t = Wild('t')
#            r = expr.match(c*t)
            term = S.One
            if expr.is_Mul:
                for i in expr.args:
                    if i.has(x):
                        term *= i
            elif expr.has(x):
                term = expr
            return term

        expr = expand_mul(expr)
        if expr.is_Add:
            for term in expr.args:
                if _remove_coefficient(term, x) in exprs:
                    pass
                else:
                    exprs.add(_remove_coefficient(term, x))
                    exprs = exprs.union(_get_trial_set(term, x, exprs))
        else:
            term = _remove_coefficient(expr, x)
            tmpset = exprs.union(set([term]))
            oldset = set([])
            while tmpset != oldset:
                # If you get stuck in this loop, then _test_term is probably broken
                oldset = tmpset.copy()
                expr = expr.diff(x)
                term = _remove_coefficient(expr, x)
                if term.is_Add:
                    tmpset = tmpset.union(_get_trial_set(term, x, tmpset))
                else:
                    tmpset.add(term)
            exprs = tmpset
        return exprs



    retdict['test'] = _test_term(expr, x)
    if retdict['test']:
        # Try to generate a list of trial solutions that will have the undetermined
        # coefficients.  Note that if any of these are not linearly independent
        # with any of the solutions to the homogeneous equation, then they will
        # need to be multiplied by sufficient x to make them so.  This function
        # DOES NOT do that (it doesn't even look at the homogeneous equation).
        retdict['trialset'] = _get_trial_set(expr, x)

    return retdict

def ode_nth_linear_constant_coeff_variation_of_parameters(eq, func, order, match):
    r"""
    Solves an nth order linear differential equation with constant
    coefficients using the method of undetermined coefficients.

    This method works on any differential equations of the form
    f(x)^(n) + a_(n-1)*f(x)^(n-1) + ... + a1*f'(x) + a0*f(x) = P(x).

    This method works by assuming that the particular solution takes the
    form Sum(c_i(x)*y_i(x), (x, 1, n)), where y_i is the ith solution to
    the homogeneous equation.  The solution is then solved using
    Wronskians and Cramer's Rule.  The particular solution is given by
    Sum(Integral(W_i(x)/W(x), x)*y_i(x), (x, 1, n)), where W(x) is the
    Wronskian of the fundamental system (the system of n linearly
    independent solutions to the homogeneous equation), and W_i(x) is
    the Wronskian of the fundamental system with the ith column replaced
    with [0, 0, ..., 0, P(x)].

    This method is general enough to solve any nth order inhomogeneous
    linear differential equation with constant coefficients, but
    sometimes SymPy cannot simplify the Wronskian well enough to
    integrate it.  If this method hangs, try using the
    'nth_linear_constant_coeff_variation_of_parameters_Integral' hint
    and simplifying the integrals manually.  Also, prefer using
    'nth_linear_constant_coeff_undetermined_coefficients' when it
    applies, because it doesn't use integration, making it faster and
    more reliable.

    Warning, using simplify=False with
    'nth_linear_constant_coeff_variation_of_parameters' in dsolve()
    may cause it to hang, because it will not attempt to simplify
    the Wronskian before integrating.  It is recommended that you only
    use simplify=False with
    'nth_linear_constant_coeff_variation_of_parameters_Integral' for
    this method, especially if the solution to the homogeneous
    equation has trigonometric functions in it.

    == Example ==
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> f = Function('f')
        >>> pprint(dsolve(f(x).diff(x, 3) - 3*f(x).diff(x, 2) +
        ... 3*f(x).diff(x) - f(x) - exp(x)*log(x), f(x),
        ... hint='nth_linear_constant_coeff_variation_of_parameters'))
               /             3 /11   log(x)\       2\  x
        f(x) = |C1 + C2*x - x *|-- - ------| + C3*x |*e
               \               \36     6   /        /

    == References ==
        - http://en.wikipedia.org/wiki/Variation_of_parameters
        - http://planetmath.org/encyclopedia/VariationOfParameters.html
        - M. Tenenbaum & H. Pollard, "Ordinary Differential Equations",
          Dover 1963, pp. 233

    """
    gensol = ode_nth_linear_constant_coeff_homogeneous(eq, func, order, match,
        returns='both')
    match.update(gensol)
    return _solve_variation_of_parameters(eq, func, order, match)

def _solve_variation_of_parameters(eq, func, order, match):
    """
    Helper function for the method of variation of parameters.

    See the ode_nth_linear_constant_coeff_undetermined_coefficients()
    docstring for more information on this method.

    match should be a dictionary that has the following keys:
    'list' - A list of solutions to the homogeneous equation, such as
         the list returned by
         ode_nth_linear_constant_coeff_homogeneous(returns='list')
    'sol' - The general solution, such as the solution returned by
        ode_nth_linear_constant_coeff_homogeneous(returns='sol')


    """
    x = func.args[0]
    f = func.func
    r = match
    psol = 0
    gensols = r['list']
    gsol = r['sol']
    wr = wronskian(gensols, x)

    if r.get('simplify', True):
        wr = simplify(wr) # We need much better simplification for some ODEs.
                          # See issue 1563, for example.

        # To reduce commonly occuring sin(x)**2 + cos(x)**2 to 1
        wr = trigsimp(wr, deep=True, recursive=True)
    if not wr:
        # The wronskian will be 0 iff the solutions are not linearly independent.
        raise NotImplementedError("Cannot find " + str(order) + \
        " solutions to the homogeneous equation nessesary to apply " + \
        "variation of parameters to " + str(eq) + "(Wronskian == 0)")
    if len(gensols) != order:
        raise NotImplementedError("Cannot find " + str(order) + \
        " solutions to the homogeneous equation nessesary to apply " + \
        "variation of parameters to " + str(eq) + "(number of terms != order)")
    negoneterm = (-1)**(order)
    for i in gensols:
        psol += negoneterm*C.Integral(wronskian(filter(lambda x: x != i, \
        gensols), x)*r['b']/wr, x)*i/r[order]
        negoneterm *= -1

    if r.get('simplify', True):
        psol = simplify(psol)
        psol = trigsimp(psol, deep=True)
    return Eq(f(x), gsol.rhs + psol)

def ode_separable(eq, func, order, match):
    r"""
    Solves separable 1st order differential equations.

    This is any differential equation that can be written as
    P(y)*dy/dx = Q(x). The solution can then just be found by
    rearranging terms and integrating:
    Integral(P(y), y) = Integral(Q(x), x). This hint uses separatevars()
    as its back end, so if a separable equation is not caught by this
    solver, it is most likely the fault of that function. seperatevars()
    is smart enough to do most expansion and factoring necessary to
    convert a separable equation F(x, y) into the proper form P(x)*Q(y).
    The general solution is::

        >>> from sympy import *
        >>> x = Symbol('x')
        >>> a, b, c, d, f = map(Function, ['a', 'b', 'c', 'd', 'f'])
        >>> pprint(dsolve(Eq(a(x)*b(f(x))*f(x).diff(x), c(x)*d(f(x))), f(x),
        ... hint='separable_Integral'))
             f(x)
           /                  /
          |                  |
          |  b(y)            | c(x)
          |  ---- dy = C1 +  | ---- dx
          |  d(y)            | a(x)
          |                  |
         /                  /

    == Example ==
    ::
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> f = Function('f')
        >>> pprint(dsolve(Eq(f(x)*f(x).diff(x) + x, 3*x*f(x)**2), f(x),
        ... hint='separable'))
            /       2   \         2
        -log\1 - 3*f (x)/        x
        ----------------- = C1 - --
                6                2

    == Reference ==
        - M. Tenenbaum & H. Pollard, "Ordinary Differential Equations",
          Dover 1963, pp. 52

    """
    x = func.args[0]
    f = func.func
    C1 = Symbol('C1')
    r = match # {'m1':m1, 'm2':m2, 'y':y}
    return Eq(C.Integral(r['m2']['coeff']*r['m2'][r['y']]/r['m1'][r['y']],
        (r['y'], None, f(x))), C.Integral(-r['m1']['coeff']*r['m1'][x]/
        r['m2'][x], x)+C1)

