# Parsing Guide

## Parsing Rules

rubi is originally written and tested on mathematica. All code and tests in mathematica format is publicly available on it's [website](http://www.apmaths.uwo.ca/~arich/).

The rubi code is parsed into sympy format in 2 steps:

#### Step 1

First the downvalues of rules are generated.
* First download the rubi folder from rubi's site and open `rubi.nb`(name may vary from version to version).

* Now run the following mathematica's code to generate the `DownValues` 

```Mathematica
LoadRules[filename_String] :=
  Module[{object},
  object=PrintTemporary["Loading "<>filename<>".m..."];
  Get[RulesDirectory<>filename<>".m"];
  NotebookDelete[object];
  Null]

inputFiles = {"9.1 Integrand simplification rules", "1.1 Linear products", "1.2 Quadratic products", "1.3 Binomial products", 
"1.4 Trinomial products", "1.5 Miscellaneous algebraic functions", "9.3 Piecewise linear functions", "2 Exponentials",
"3 Logarithms", "4.1 Sine", "4.2 Tangent", "4.3 Secant", "4.4 Miscellaneous trig functions", "5 Inverse trig functions",
"6 Hyperbolic functions", "7 Inverse hyperbolic functions", "8 Special functions", "9.4 Miscellaneous integration rules"}

outputFiles = {"Integrand_simplification.txt", "Linear_products.txt", "Quadratic_products.txt", "Binomial_products.txt",
"Trinomial_products.txt", "Miscellaneous_algebra.txt", "Piecewise_linear.txt", "Exponentials.txt", "Logarithms.txt",
"Sine.txt", "Tangent.txt", "Secant.txt", "Miscellaneous_trig.txt", "Inverse_trig.txt", "Hyperbolic.txt",
"Inverse_hyperbolic.txt", "Special_functions.txt", "Miscellaneous_integration.txt"}

ShowSteps = False
LoadRules["ShowStep routines"];
LoadRules["Integration utility functions"];

For[i = 1, i < 19, i++, Clear[Int]; Int::usage="Int [expn, var]"; LoadRules[inputFiles[[i]]];
Unprotect[Sinc]; Sinc[u_] := Sin[u]/u; Protect[Sinc];
FixIntRules[]; Export[outputFiles[[i]], ToString@FullForm@DownValues@Int];
]
```

(Note that you can change the items in inputFiles according to the files present in rubi).
Now there are new files as in `outputFiles` containing the `DownValues`.

#### Step 2

Second step is to generate sympy code from `DownValues`. For this there is a function written in `generate_rules.py`.
Run the following code in python terminal:

```python
>>> from sympy.integrals.rubi.parsetools.generate_rules import generate_rules_from_downvalues
>>> generate_rules_from_downvalues()
```

You can find the parsed rules. Also a file `constraints.py` containing all constraints.

(Note : Be careful with the name of files. The output files name from step1 is exactly same as input files name in step2)

-------------------------------------------------

## Parsing Tests

rubi contains a large set of test cases in mathematica format. Those need to be parsed in sympy format for testing. This is done in two steps.

#### Step1

In first step, we need to get the `FullForm` of tests.

* Download test files from official rubi [website](http://www.apmaths.uwo.ca/~arich/IntegrationProblems/MathematicaSyntaxFiles/MathematicaSyntaxFiles.html) in mathematica format.

* Open a mathematica notebook and run the following mathematica code.(Here, we are assuming the name of downloaded test file to be `testMath.m`. This needs to be changed as per the situation.)

```Mathematica
stream = OpenWrite["test_1.m"];
WriteString[stream, Import["testMath.m", "HeldExpressions"] // FullForm];
Close[stream1];
```

The above code writes `FullForm` of test cases in `test_1.m`

#### Step 2

Next, we need to parse `test_1.m` in sympy format.
Open a python terminal and run the following code: 

```python
>>> from sympy.integrals.rubi.parsetools.generate_tests import generate_test_file
>>> generate_test_file()
```

The above code writes tests in sympy format in `parsed_tests.py`. File names in `generate_tests` should be changed as per the situation.

Note: Current test suite in sympy is not all parsed through these above steps. But it works well and is tested for `special_functions`. `test_error_functions` in `test_special_functions.py` has been parsed through the above steps.

-----------------------

### References

* [1] http://reference.wolfram.com/language/ref/FullForm.html
* [2] http://reference.wolfram.com/language/ref/DownValues.html
