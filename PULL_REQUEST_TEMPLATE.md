<!--
Thank you so much for your PR!  To help us review your contribution, please
consider the following points:

- Your title above should be a short description of what was changed. Do not
  include the issue number in the title.

- Development information is available at
  https://github.com/sympy/sympy/wiki/Development-workflow.

- Do not create the PR out of master, but out of a separate branch.

- Try to keep the number of commits to a minimum. Use rebase and
  force pushing in git to avoid adding unnecessary commits.

- If this is your first commit, see:
  https://github.com/sympy/sympy/wiki/Development-workflow#add-your-name-and-email-address-to-the-mailmap-file

- If you are contributing fixes to docstrings, please pay attention to
  https://docs.sympy.org/dev/guides/contributing/documentation-style-guide.html.
  In particular, note the difference between using single backquotes,
  double backquotes, and asterisks in the markup.

- Please include tests of any new feature or bug fix. The code in the tests
  should not execute correctly in a version of SymPy without your code.

Note that all reviewing is done an a volunteer basis, so it can sometimes take
time before you get feedback. If you have not heard back in two weeks, feel
free to ping @sympy/developers-with-push-access-to-everything
Do not ping individual persons to get a review, nor request reviews.
Also note that any reviewer will be notified of changes in your PR, so there
is no need to ping after a change.
-->

#### References to other Issues or PRs
<!-- If this pull request fixes an issue, write "Fixes #NNNN" in that exact
format, e.g. "Fixes #1234" (see https://tinyurl.com/auto-closing for more
information). Also, please write a comment on that issue linking back to this
pull request once it is open. -->


#### Brief description of what is fixed or changed


#### Other comments


#### Release Notes

<!-- Write the release notes for this release below between the BEGIN and END
statements. The basic format is a bulleted list with the name of the subpackage
and the release note for this PR. For example:

* solvers
  * Added a new solver for logarithmic equations.

* functions
  * Fixed a bug with log of integers. Formerly, `log(-x)` incorrectly gave `-log(x)`.

* physics.units
  * Corrected a semantical error in the conversion between volt and statvolt which
    reported the volt as being larger than the statvolt.

or if no release note(s) should be included use:

NO ENTRY

See https://github.com/sympy/sympy/wiki/Writing-Release-Notes for more
information on how to write release notes. The bot will check your release
notes automatically to see if they are formatted correctly. -->

<!-- BEGIN RELEASE NOTES -->

<!-- END RELEASE NOTES -->
