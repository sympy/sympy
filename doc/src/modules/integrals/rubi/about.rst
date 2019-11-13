About Rule Based Integrator
===========================

About it
--------

A RUle Based Integrator nicknamed RUBI is a module in sympy.integrals. It is an
implementation of more than 6000 rules to determine antiderivative of a wide
variety of mathematical expressions. RUBI utilizes a set of well defined rules
which makes it smart to present the results in a more symmetric and simplified
manner. The Pattern Matcher (MatchPy) along with transformation rules play an
important role in making RUBI more robust and fast. As it works on the given
set of rules. MatchPy is a pattern matching library in Python which has matching
capabilities similar to Mathematica. But because MatchPy is only implemented in
Python3.6 so RUBI is not supported by Python versions less than 3.6.

Note
----

Please note that this module is experimental. Therefore, API may change unexpectedly
between SymPy versions.

References
----------
http://www.apmaths.uwo.ca/~arich/
