Pattern Matcher
===============

A general purpose repository might require thousands or even tens of thousands 
of rules. Obviously sequentially searching a list of that many rules to find a 
match would be unacceptably slow. Even having a separate list of rules for each
built-in function or operator is insufficient, since some functions may have a 
large number of rules associated with it (e.g. our integrator requires over a 
1000 rules). 

Matchpy
-------

**MatchPy** implements pattern matching in python. It is similar to the 
implementation in Mathematica. A detailed example of how you can use matchpy 
can be found in its documentation (see the link below). In addition to the 
basic matching algorithm, there are data structures that can be used for more
efficient many-to-one matching like the ManyToOneMatcher and the 
DiscriminationNet.

References
----------

https://matchpy.readthedocs.io/en/latest/