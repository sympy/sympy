## Fixed mishandling of power function in the parsing.mathematica module
<!-- Your title above should be a short description of what
was changed. Do not include the issue number in the title. -->

#### References to other Issues or PRs
Fixes #24150
<!-- If this pull request fixes an issue, write "Fixes #NNNN" in that exact
format, e.g. "Fixes #1234" (see
https://tinyurl.com/auto-closing for more information). Also, please
write a comment on that issue linking back to this pull request once it is
open. -->


#### Brief description of what is fixed or changed

Fixed wrong behavior of function **parse_mathematica** in module ***parsing.mathematica***.
In function parse changed statement ***_mathematica_op_precedence*** as follow:

*Rule of power* transferred **before** *rule of prefix operators*.
```python
#rule of power transferred before rule of prefix operators
(INFIX, RIGHT, {"^": "Power"}),
(PREFIX, None, {"-": lambda x: MathematicaParser._get_neg(x),
                "+": lambda x: x}),
```

#### Other comments


#### Release Notes

<!-- Write the release notes for this release below between the BEGIN and END
statements. The basic format is a bulleted list with the name of the subpackage
and the release note for this PR. For example:

* solvers
  * Added a new solver for logarithmic equations.

* functions
  * Fixed a bug with log of integers.

or if no release note(s) should be included use:

NO ENTRY

See https://github.com/sympy/sympy/wiki/Writing-Release-Notes for more
information on how to write release notes. The bot will check your release
notes automatically to see if they are formatted correctly. -->

<!-- BEGIN RELEASE NOTES -->

* parsing
  * functions
    * Fixed a wrong behavior (parse of power negative arguments) of function parse_mathematica.
    * Make `parse_mathematica("2^-10*x")` give `x/1024` instead of `2**(-10*x)`.

<!-- END RELEASE NOTES -->
