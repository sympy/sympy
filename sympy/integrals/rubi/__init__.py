'''
Rule Based Integration(RUBI) module in sympy uses set of transformation
rules to integrate an expression. All the transformation rules are compiled as a
discrimination-net which helps in matching expression with the rule efficiently.

Due to large number of rules, the module would normally take lot of time to load.
Hence, it is better to use Rubi while doing multipe integrations.

Note: This module has dependency on MatchPy library.

Example
=======
>>> from sympy.integrals.rubi.rubi import rubi_integrate
>>> from sympy.abc import x
>>> rubi_integrate(1/x, x)
log(x)

It is possible to know which rules are matching with the expression by using
`get_matching_rule_definition()` function

>>> from sympy.integrals.rubi.rubi import get_matching_rule_definition
>>> get_matching_rule_definition(1/x, x)
Rule matching:
/Users/parsoyaarihant/sympy/sympy/integrals/rubi/rules/linear_products.py
On line:  21
    rule1 = ReplacementRule(pattern1, lambda x : Log(x))
<BLANKLINE>
Pattern matching:
(Pattern(Integral(u*(a + b*x**n)**p, x), constraints=(CustomConstraint(FreeQ(a, x)), CustomConstraint(FreeQ(a, x)), CustomConstraint(FreeQ(a, x)), CustomConstraint(FreeQ(a, x)), CustomConstraint(UNKNOWN))), <function integrand_simplification.<locals>.<lambda> at 0x10574b7b8>, [0, 1, 2, 3, 4])
{x ↦ x}

TODO
====
* Use code generation to implement all rules.
* Testing of all the tests from rubi test suit. See: http://www.apmaths.uwo.ca/~arich/IntegrationProblems/MathematicaSyntaxFiles/MathematicaSyntaxFiles.html
* Add support for `Piecewise` functions.

References
==========
[1] http://www.apmaths.uwo.ca/~arich/
[2] https://github.com/sympy/sympy/issues/7749
[3] https://github.com/sympy/sympy/pull/12978
'''
