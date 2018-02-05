'''
Rule Based Integration(RUBI) module in sympy uses set of transformation
rules to integrate an expression. All the transformation rules are compiled as a
discrimination-net which helps in matching expression with the rule efficiently.

Due to large number of rules, the module would normally take lot of time to load.
Hence, it is better to use Rubi while doing multiple integrations.

Rules are taken from Rubi version 4.10.8.

Note: This module has dependency on MatchPy library.

TODO
====
* Use code generation to implement all rules.
* Testing of all the tests from rubi test suit. See: http://www.apmaths.uwo.ca/~arich/IntegrationProblems/MathematicaSyntaxFiles/MathematicaSyntaxFiles.html
* Add support for `Piecewise` functions.

Debugging
=========
When an integration is not successful. We can see which rule is matching the
expression by using `get_matching_rule_definition()` function. We can cross-check
if correct rule is being applied by evaluating the same expression in Mathematica.
If the applied rule is same, then we need to check the `ReplacementRule` and
the utility functions used in the `ReplacementRule`.

Parsing Rules
=============
I have included the parser in /parsetools folder. The parser takes
`FullForm[DownValues[]]` of Mathematica rules as input.

References
==========
[1] http://www.apmaths.uwo.ca/~arich/
[2] https://github.com/sympy/sympy/issues/7749
[3] https://github.com/sympy/sympy/pull/12978
'''
