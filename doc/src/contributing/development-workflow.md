# Development workflow

## Introduction

In [SymPy] we encourage collaborative work.

Everyone is welcome to join and to implement new features, fix bugs, give
general advice, etc. Also, we try to discuss everything and to review each
other's work so that many eyes can see more thus raising the quality.

General discussion takes place on the [sympy@googlegroups.com] mailing list and
in the [issue tracker]. Some discussion also takes place on Gitter (our channel is
[Gitter sympy/sympy]).

As some of you already know, software development is not just coding. Many
non-coding tasks have to be done in order to produce *good* code. For
example: setting up infrastructure, designing, testing, documenting,
assisting new developers (we are doing it here), and of course programming.

As already said above, we review changes. In short, each change is first
reviewed by other developers and only when it is approved is the code pushed
in. This is important to produce the highest quality code, and to ensure the
functions in SymPy are correct.

And just as it takes effort to write good and clear code, reviewing other's
work needs effort too. There are good practices how to do this so that
reviewing is fun for both the author and the reviewer. We try to follow these
good practices, and we'll try to show you how to follow them too.

When reviewing other's patches you *learn* a lot, so why not participate
as a reviewer too? Anyone regardless of technical skill can help review code,
and it's an excellent way for newcomers to learn about SymPy's development
process and community.

## How to contribute to the code base

License: [New BSD License] (see the [LICENSE] file for details) covers all files in the SymPy repository unless stated otherwise.
By submitting your patch, you are agreeing to let your work be licensed under this.

There are a few ways to create and send a patch.

The best way is to send a GitHub pull request against the [sympy/sympy] repository. We'll review it and push it in.
The GitHub pull request is the preferred method, because it makes it easy for us to review and push the code in.

The basic workflow is as follows:

1. [Create your environment](#create-your-environment), if it was not created earlier.
2. Pick an issue to fix.
3. Create a new branch.
4. Modify code and/or create tests of it.
5. Be sure that all tests of SymPy pass.
6. Only then commit the changes.
7. Create a patch file or pull request for GitHub.

After you've submitted your patch, it will be reviewed by reviewers.
They may request you to make changes, in which case you may need to:

- Update your pull request
- Synchronize with master sympy/sympy

All these are described in [Workflow process](#workflow-process), but before
you read that, it would be useful to acquaint yourself with [Coding
conventions in SymPy](#coding-conventions-in-sympy).

If you have any questions you can ask them on the [mailing list].

## Coding conventions in SymPy

### Standard Python coding conventions

Follow the standard Style Guide for Python Code when writing code for SymPy, as explained at the following URLs:

- <https://www.python.org/dev/peps/pep-0008>
- <https://www.python.org/dev/peps/pep-0257>

In particular,

- Use 4 spaces for indentation levels.

- Use all lowercase function names with words separated by
  underscores. For example, you are encouraged to write Python
  functions using the naming convention

  ```
  def set_some_value()
  ```

  instead of the CamelCase convention.

- Use CamelCase for class names and major functions that create
  objects, e.g.

  ```
  class PolynomialRing(object)
  ```

Note, however, that some functions do have uppercase letters where it makes
sense. For example, for matrices they are LUdecomposition or T (transposition)
methods. Also, classes that represent mathematical functions use the case
convention most commonly used for that function. For example, functions like
`sin` and `exp` are implemented as classes, but they are spelled in all
lowercase.

### Documentation strings

For information on writing and formatting docstrings, please see the [SymPy Documentation Style Guide](documentation-style-guide).

### Python 3

SymPy only supports Python 3. If you are using Python 2, SymPy will not work, and you will need to install Python 3 first in order to use it.

## Create your environment

The first step to contributing to the code base is creating your development environment.
This is done only once.

### Install Anaconda

The easiest way to install Python is to install Anaconda. <https://www.anaconda.com/products/distribution>. This will also install most of the dependencies you need to use SymPy.

### Install mpmath

SymPy has a hard dependency on the [mpmath](http://mpmath.org/) library
(version >= 0.19). You should install it first. If you installed Anaconda, it
already comes with mpmath, so you do not need to worry about this step.
Otherwise you can install it with conda

```
conda install mpmath
```

or pip

```
pip install mpmath
```

### Optional Dependencies

Depending on which part of SymPy you plan to contribute to, you may need to
install one of SymPy's optional dependencies. The optional dependencies used
by SymPy are listed on [the dependencies page](dependencies.md). Most code
contributions will not require installing these dependencies. You should not
worry about installing an optional dependency unless you need it, as some of
SymPy's optional dependencies can be hard to install and are only required for
a small subset of the library.

If you are contributing documentation, you may wish to build the HTML
documentation, in which case, you will need to [install the dependencies to
build the documentation](build-docs.rst).

### Set up git

#### Install git

**Linux-like systems**:

Install git via your native package management system:

```
yum install git
```

or:

```
sudo apt-get install git
```

**Windows and macOS**:

The easiest way to get git is to download [GitHub
desktop](https://desktop.github.com/), which will install git, and also
provide a nice GUI (this tutorial will be based on the command line
interface). Note, you may need to go into the GitHub preferences and choose
the "Install Command Line Tools" option to get git installed into the
terminal.

If you do decide to use the GitHub GUI, you should make sure that any "sync
does rebase" option is disabled in the settings.

#### Configure git settings

Git tracks who makes each commit by checking the user’s name and email.
In addition, we use this info to associate your commits with your GitHub account.

To set these, enter the code below, replacing the name and email with your own (`--global` is optional).:

```
git config --global user.name "Firstname Lastname"
git config --global user.email "your_email@youremail.com"
```

The name should be your actual name, not your GitHub username.

These global options (i.e. applying to all repositories) are placed in `~/.gitconfig`.
If you want, you can edit this file to add setup colors and some handy shortcuts:

```
[user]
    name = Firstname Lastname
    email = your_email@youremail.com

[color]
    diff  = auto
    status= auto
    branch= auto
    interactive = true

[alias]
    ci = commit
    di = diff --color-words
    st = status
    co = checkout
    log1 = log --pretty=oneline --abbrev-commit
    logs = log --stat
```

#### Tune bash prompt

It can be convenient in future to tune the bash prompt to display the current git branch.

The easiest way to do it, is to add the snippet below to your .bashrc or .bash_profile:

```
PS1="[\u@\h \W\$(git branch 2> /dev/null | grep -e '\* ' | sed 's/^..\(.*\)/{\1}/')]\$ "
```

But better is to use `git-completion` from the `git` source. This also has the advantage of adding tab completion to just about every git command. It also includes many other useful features, for example,
promptings. To use `git-completion`, first download the `git` source code (about 27 MiB), then copy
the file to your profile directory:

```
git clone git://git.kernel.org/pub/scm/git/git.git
cp git/contrib/completion/git-completion.bash ~/.git-completion.sh
```

Read instructions in '~/.git-completion.sh'

Note that if you install git from the package manager in many Linux distros, this file is already installed for you.  You can check if it is installed by seeing if tab completion works on git commands (try, e.g., `git commi<TAB>`, or `git log --st<TAB>`). You can also check if the PS1 commands work by doing something like:

```
PS1='\W $(__git_ps1 "%s")\$ '
```

And your command prompt should change to something like:

```
sympy master$
```

Note, it is important to define your PS1 using single quotes ('), not double quotes ("), or else bash will not update the branch name.

#### Create GitHub account

As you are going to use [GitHub]  you should have a GitHub account. If you have not one yet then sign up at:

- <https://github.com/signup/free>

#### Set up SSH keys

To establish a secure connection between your computer and GitHub see detailed instructions in <https://help.github.com/en/articles/set-up-git> or at <https://help.github.com/en/articles/adding-a-new-ssh-key-to-your-github-account>.

If you have any problems with SSH access to GitHub, read the troubleshooting
instructions at <https://help.github.com/en/articles/troubleshooting-ssh>, or
ask us on the [mailing list].

#### Fork SymPy project

Create your own *fork* of the SymPy project (if you have not yet). Go to the SymPy GitHub repository:

- <https://github.com/sympy/sympy>

and click the “Fork” button.

Now you have your own repository for the SymPy project. The address of the forked project will look something like:

- <https://github.com/<your-github-username>/sympy>

where `<your-github-username>` is your GitHub username.

#### Clone SymPy

On your machine browse to where you would like to store SymPy, and clone (download) the latest
code from SymPy's original repository (about 77 MiB):

```
git clone git://github.com/sympy/sympy.git
cd sympy
```

Then assign your read-and-write repo to a remote called "github" (replace
`<your-github-username>` with your GitHub username):

```
git remote add github git@github.com:<your-github-username>/sympy.git
```

For more information about GitHub forking and tuning see:
<https://help.github.com/en/articles/about-pull-requests>, <https://help.github.com/en/articles/fork-a-repo>, and <https://help.github.com/en/articles/set-up-git>


### Set up virtual environments

You may want to take advantage of using virtual environments to isolate your development version of SymPy from any system wide installed versions, e.g. from `apt-get install python-sympy`. There are two leading virtual environment tools, [virtualenv](https://virtualenv.pypa.io) and [conda](http://conda.pydata.org/). Conda comes with Anaconda, which is what we recommended to install above. Here is an example of using conda to create a virtual environment:

```
conda create -n sympy-dev python=3 mpmath flake8
```

Note that you should *not* include `sympy` as one of the packages you install
into this environment. That would install the released version of SymPy, but
you want to work with the development version of SymPy.

You now have a environment that you can use for testing your development copy
of SymPy.

Now activate the environment:

```
conda activate sympy-dev
```

Then, from the fork you can run the SymPy tests (these may take about an hour
to run depending on your machine):

```
(sympy-dev)$ ./bin/test
```

You may also want to try out the Flake8 style enforcement tool. If everything works fine, there should be no output from `flake8`. This command may take a few minutes to complete:

```
(sympy-dev)$ flake8 sympy
```

Optionally, you can install SymPy into the environment if you wish. This step
is not required as you can always test your changes by running Python from
within the `sympy` directory. This is only required if you wish to use the
SymPy development version from any location on your filesystem.

```
(sympy-dev)$ pip install -e .
```

If you prefer virtualenv instead of conda, the process is similar.

## Workflow process

### Pick an issue to fix

The best way to start with the main code base is to fix some existing bugs. Peruse the [Easy to fix issues] in the issue tracker and see if one interests you. If you'd like to try to fix it, then create a message in the issue saying that you'd like to work on it. If it isn't clear how to fix it, ask for suggestions on how to do it in the issue itself, on the [mailing list], or on [Gitter sympy/sympy].

SymPy's code is organized into Python packages and modules. The core code is
in the `sympy/core` directory and other packages in the sympy directory have
more specific code. For example, `sympy/printing` contains the code that
handles how SymPy objects are printed to the terminal and Jupyter.

If you are making a change that does not yet have an issue, it is not required
to open an issue first. This is only necessary if you feel you want to discuss
the change before you make a pull request, for example, if you're not sure if
something is actually a bug, or if you're not sure if a new feature is in
scope. It's also fine to just open a pull request and discuss it there.
Discussions are easier to have when there is actual code to talk about, so a
pull request is preferable if you have changes, even if they aren't fully
ready to merge yet.

### Create a new branch

The first thing to do before making changes to the code is to make a branch in git.

Remember, **never begin your work in master**. `master` should only be used to
pull in upstream changes from the main sympy/sympy repo. If you commit to
`master`, it will be difficult to pull these changes in, and it will also be
difficult if you wish to make more than one pull request at a time.

To create and checkout (that is, make it the working branch) a new branch, say
`fix-solve-bug`, run

```
git branch fix-solve-bug
git checkout fix-solve-bug
```

or in one command using

```
git checkout -b fix-solve-bug
```

To view all branches, with your current branch highlighted, type:

```
git branch
```

And remember, **never type the following commands in master**: `git merge`, `git add`, `git commit`, `git rebase`.
If you had made some commits to your local master by accident, you would have to hard reset to drop the commits.

#### Branch names

Use a short, easy to type branch name that somehow relates to the changes,
e.g., `fix-solve-bug`. Avoid using issue numbers in branch names as these are
not easy to type (most SymPy issue numbers are 5 digits long) and they won't
really be indicative of what the change is about without looking up the issue.

### Modify code

Do not forget that all new functionality should be tested, and all new methods, functions, and classes should have [doctests] showing how to use them.

Keep in mind, doctests are *not* tests. Think of them as examples that happen to be tested. Some key differences:

- write doctests to be informative; write regular tests to check for regressions and corner cases.
- doctests can be changed at any time; regular tests should not be changed.

In particular, we should be able to change or delete any doctest at any time if it makes the docstring better to understand.

### Be sure that all tests of SymPy pass

To ensure everything stays in shape, let’s see if all tests pass:

```
./bin/test
./bin/doctest
```

If you are on Windows, you may need to run the commands like `python bin/test` instead.

Each command will show a *DO NOT COMMIT* message if any of the tests it runs does not pass.

You can also run `bin/test --slow`, to run the slow tests. These take a long
time to run.

Code quality (unwanted spaces and indents) are checked by *./bin/test* utilities too. But you can separately run this test with the help of this command:

```
./bin/test quality
```

If you have trailing whitespace it will show errors. This one will fix unwanted spaces:

```
./bin/strip_whitespace <file>
```

It's best to configure your editor settings to automatically remove trailing
whitespace when you save so that you never have to worry about this.

If you want to test only one set of tests try:

```
./bin/test sympy/concrete/tests/test_products.py
```

If you want to test only one specific test try:

```
./bin/test sympy/concrete/tests/test_products.py -k test_one
```

But remember that all tests should pass before committing.

This includes the style checks done by the Flake8 tool. You can run these additional tests using the following command:

```
flake8 sympy
```

Just as with SymPy's own test suite, it's also possible to restrict `flake8` checking to a single file by appending its name on the command line.

Note that all tests will be run when you make your pull request automatically
by GitHub Actions, so do not worry too much about running every possible test.
You can usually just run:

```
./bin/test mod
./bin/doctest mod
```

where `mod` is the name of the module that you modified.

### Commit the changes

You can check what files are changed:

```
git status
```

Check total changes:

```
git diff
```

If you created any new files, add them with:

```
git add new_file.py
```

You are ready to commit changes locally. A commit also contains a `commit
message` which describes it.  See the next section for guidelines on writing
good commit messages. Type:

```
git commit
```

An editor window will appear automatically in this case. In Linux, this is vim by default. You
can change what editor pops up by changing the `$EDITOR` shell variable.

Also with the help of option `-a` you can tell the command `commit` to automatically stage files
that have been modified and deleted, but new files you have not told git about will not be
affected, e.g.,:

```
git commit -a
```

If you want to stage only part of your changes, you can use the interactive commit feature.  Just type:

```
git commit --interactive
```

and choose the changes you want in the resulting interface.

#### Deleting junk files

A lot of editors can create some configuration files, binary files, or temporary files
in your SymPy directory, which should be removed before merging your commits.

Tracking down individual files can be cumbersome.

You may think of using `.gitignore`, however, editing the `.gitignore` itself
would have the agreement from the community.

Using `.git/info/exclude` would be the best, because it is only applied locally.

<https://stackoverflow.com/questions/22906851/when-would-you-use-git-info-exclude-instead-of-gitignore-to-exclude-files>

<https://help.github.com/en/articles/ignoring-files>

(development-workflow-commit-messages)=
#### Writing commit messages

The commit message has two parts: a title (first line) and the body. The two
are separated by a blank line.

Commit messages summarize what the commit does. Just as with the code, your
commit messages will become a permanent part of the project git history. So
you should put some effort into making them high quality. Commit messages are
intended for human readers, both for people who will be reviewing your code
right now, and for people who might come across your commit in the future
while researching some change in the code. Thus, include information that
helps others understand your commit here, if necessary.

Tools like `git shortlog` or even GitHub only show the first line of the commit by
default, so it is important to convey the most important aspects of the commit in the first line.

- Keep the first line 71 characters or less and subsequent lines to 78
  characters or less. This allows the one-line form of the log to display the
  summary without wrapping.

- **Make sure to leave a blank line after the summary**

- Do not end the first line with a period (full stop). Subsequent lines should
  use periods.

- Provide context for the commit if possible,

  e.g. `integrals: Improved speed of heurisch()`
  instead of just `Improved speed of heurisch()`

- Reference any relevant issue numbers. You do not need to reference the pull
  request for the change itself, but issues that are fixed should be
  referenced either by `#12345` or
  `https://github.com/sympy/sympy/issues/12345`. You should also provide a
  brief summary of an issue rather than just referring to the issue number so
  that people don't have to look around for context.

- A commit won't always be seen in the context of your branch, so it is often
  helpful to give each commit some context. This is not required, though, as
  it is not hard to look at the commit metadata to see what files were
  modified or at the commit history to see the nearby related commits.

- Use plain English. Write in complete sentences.

- Describe what actually changed. Don't just write something like `Modified
  solvers.py`. People can already see what files are modified from the commit
  diff. What the message is there for is to tell them what the diff actually
  does, so they don't have to try to figure it out.

- Try to avoid short commit messages, like "Fix", and commit messages that
  give no context, like "Found the bug". When in doubt, a longer commit
  message is probably better than a short one. Avoid using the `-m` switch to
  `git commit` to write a commit message on the command line. Rather, let it
  open your editor so you can write a longer commit message.

- Give an overview of what the commit does if it is difficult to figure out just
  from looking at the diff.

- Include other relevant information, e.g.

  - Known issues
  - A concrete example (for commits that add new features/improve performance etc.)

- Use bullet lists when suitable

- Feel free to use Unicode characters, such as output from the SymPy Unicode pretty printer.

- Use plain English

##### Example of a good commit message

Here is an example commit message from the commit
[bf0e81e12a2f75711c30f0788daf4e58f72b2a41](https://github.com/sympy/sympy/commit/bf0e81e12a2f75711c30f0788daf4e58f72b2a41),
which is part of the SymPy history:

```
integrals: Improved speed of heurisch() and revised tests

Improved speed of anti-derivative candidate expansion and solution
phases using explicit domains and solve_lin_sys(). The upside of
this change is that large integrals (those that generate lots of
monomials) are now computed *much* faster. The downside is that
integrals involving Derivative() don't work anymore. I'm not sure
if they really used to work properly or it was just a coincidence
and/or bad implementation. This needs further investigation.

Example:

In [1]: from sympy.integrals.heurisch import heurisch

In [2]: f = (1 + x + x*exp(x))*(x + log(x) + exp(x) - 1)/(x + log(x) + exp(x))**2/x

In [3]: %time ratsimp(heurisch(f, x))
CPU times: user 7.27 s, sys: 0.04 s, total: 7.31 s
Wall time: 7.32 s
Out[3]:
   ⎛ 2        x                 2⋅x      x             2   ⎞
log⎝x  + 2⋅x⋅ℯ  + 2⋅x⋅log(x) + ℯ    + 2⋅ℯ ⋅log(x) + log (x)⎠          1
──────────────────────────────────────────────────────────── + ───────────────
                             2                                      x
                                                               x + ℯ  + log(x)

Previously it took 450 seconds and 4 GB of RAM to compute.
```

##### Co-Author

Occasionally, there can be multiple people working as a team for one PR,
or you have applied some suggestions from the community.

For these cases, you may use co-author feature of GitHub by adding

```
Co-authored-by: NAME NAME@EXAMPLE.COM
Co-authored-by: AUTHOR-NAME ANOTHER-NAME@EXAMPLE.COM
```

to the bottom of the commit message. See
https://help.github.com/en/articles/creating-a-commit-with-multiple-authors.

### Create a patch file or pull request for GitHub

Be sure that you are in your own branch, and run:

```
git push github fix-solve-bug
```

This will send your local changes to your fork of the SymPy repository.

When you push for the first time, it will print a link that you can use to
open a pull request, like

```
Enumerating objects: 32, done.
Counting objects: 100% (32/32), done.
Delta compression using up to 12 threads
Compressing objects: 100% (27/27), done.
Writing objects: 100% (32/32), 6.42 KiB | 2.14 MiB/s, done.
Total 32 (delta 8), reused 22 (delta 5), pack-reused 0
remote: Resolving deltas: 100% (8/8), done.
remote:
remote: Create a pull request for 'fix-solve-bug' on GitHub by visiting:
remote:      https://github.com/<your-github-username>/sympy/pull/new/fix-solve-bug
remote:
To https://github.com/<your-github-username>/sympy.git
 * [new branch]      fix-solve-bug -> fix-solve-bug
 ```

You can click the
`https://github.com/<your-github-username>/sympy/pull/new/fix-solve-bug` link
in the terminal.

Or, go to https://github.com/sympy/sympy/compare, click "compare across
forks", and select your fork and your branch in the right hand side (the arrow
should point to `base repository: sympy/sympy`, `base:master`).

After clicking `Create Pull Request`, you will be presented with a preview
page containing a textbox for the title, a textbox for the description, and
the commits that are included

Write a descriptive title for the pull request, and fill in the fields in the
description. See below for more details on what to put in the pull request
description.

You can double check that you are committing the right changes by switching to
the `Commits` tab to see which commits are included (sometimes unintended
commits can be caught this way) and switching to the `Files Changed` tab to
review the diff of all changes

When you are ready, press the `Send pull request` button. The pull request is
sent immediately and you’re taken to the main pull request discussion and
review page. Additionally, all repository collaborators and followers will see
an event in their dashboard.

Now, people will review the pull request, and you should update it with any
changes based on review comments. See [](development-workflow-update-your-pull-request).

#### Writing pull request title and description

You might feel that all your documentation work is done if you have made good commit messages.
But a good title and description will help in the review process.

- **Descriptive titles are very important.** The pull request title should
  indicate what is fixed. Pull requests with undescriptive titles will tend to
  be ignored.

  Examples of bad pull request titles are

  - "Modified solvers\.py"
  - "Fix issue #12345"

  These do indicate to reviewers what is actually changed in the pull request,
  and hence, they are likely to just pass it over instead of reviewing it.
  An example of a better pull request title is

  - "Fix a bug with solve() on transcendental functions"

- Do not put issue numbers or file names in the pull request title. Issue
  numbers should go in the description.

- Include the prefix "\[WIP\]" if you aren't ready to have the pull request
  merged and remove the prefix when you are ready

The description is a good place to:

- Show what you have done, perhaps comparing output from master with the output after your changes

- Refer to the issue that was addressed like "#1234"; that format will
  automatically create a link to the corresponding issue or pull request, e.g.
  "This is similar to the problem in issue #1234...". This format also works
  in the discussion section of the pull request.

- Use phrases like "closes #1234" or "fixed #1234" (or similar that [follow
  the auto-close
  syntax](https://help.github.com/articles/closing-issues-via-commit-messages)
  and are also [discussed
  here](https://github.com/blog/1506-closing-issues-via-pull-requests)). Then
  those other issues or pull requests will be closed automatically when your
  pull request is merged. Note: this syntax does not work in the discussion of
  the pull request. See this [quick
  guide](https://github.com/sympy/sympy/wiki/Issue-PR-Autoclosing-syntax) for
  the valid and invalid syntax for automatically closing issues from pull
  requests.

- the pull request needs a release notes entry. See
  <https://github.com/sympy/sympy/wiki/Writing-Release-Notes> on how to write
  release notes in the pull request description. The SymPy Bot will check that
  your PR has release notes automatically.

It is best to just fill out the pull request template (the text that is there
when you open a pull request). If you fill out all the sections in the
template, you will have a good pull request description.

See also [GitHub's guidelines for pull
requests](https://github.com/blog/1943-how-to-write-the-perfect-pull-request)

(development-workflow-update-your-pull-request)=
### Update your pull request

If you need to make changes to a pull request there is no need to close it.
The best way to make a change is to add a new commit in your local repository
and repeat the push command:

```
git commit
git push github fix-solve-bug
```

Note that if you do any rebasing or in any way edit your commit history, you will have to add
the `-f` (force) option to the push command for it to work:

```
git push -f github
```

You don't need to do this if you merge, which is the recommended way.

### Synchronize with master `sympy/sympy`

Sometimes, you may need to merge your branch with the upstream master. Usually
you don't need to do this, but you may need to if

- Someone tells you that your branch needs to be merged because there are
  merge conflicts.
- GitHub tells you that your branch could not be merged.
- You need some change from master that was made after you started your branch.

Note, that after cloning a repository, it has a default remote called `origin`
that points to the `sympy/sympy` repository.  And your fork remote named as
`github`. You can observe the remotes names with the help of this command:

```
$ git remote -v
github  git@github.com:mynick/sympy.git (fetch)
github  git@github.com:mynick/sympy.git (push)
origin  git://github.com/sympy/sympy.git (fetch)
origin  git://github.com/sympy/sympy.git (push)
```

As an example, consider that we have these commits in the master branch of
local git repository:

```
A---B---C        master
```

Then we have divergent branch `fix-solve-bug`:

```
A---B---C           master
         \
          a---b     fix-solve-bug
```

In the meantime the remote `sympy/sympy` master repository was updated too:

```
A---B---C---D       origin/master
A---B---C           master
         \
          a---b     fix-solve-bug
```

There are basically two ways to get up to date with a changed master: merging
and rebasing. **In general, rebasing should only be used before you make a pull
request to SymPy. Once the pull request is made your code is considered public
and the history should not be changed through rebasing. See below for the
reasons for this.**

#### Merging

Merging creates a special commit, called a "merge commit", that joins your
branch and master together:

```
A---B---C------D       origin/master
         \      \
          \      M     merge
           \    /
            a--b       fix-solve-bug
```

Note that the commits `A`, `B`, `C`, and `D` from master and the
commits `a` and `b` from `fix-solve-bug` remain unchanged. Only the new
commit, `M`, is added to `fix-solve-bug`, which merges in the new commit
branch from master.

#### Rebasing

Rebasing essentially takes the commits from `fix-solve-bug` and reapplies
them on the latest master, so that it is as if you had made them from the
latest version of that branch instead. Since these commits have a different
history, they are different (they will have different SHA1 hashes, and will
often have different content):

```
A---B---C---D---a'---b' origin/master
```

Rebasing is required if you want to edit your commit history (e.g., squash
commits, edit commit messages, remove unnecessary commits).

**In general, if you want to rebase you should only rebase before you submit a pull request.**
And this applies more for working as a team,
because rebasing breaks the commit log which other people have forked on.

However, this may not apply much if you are not working for the PR in a team.
Though, people from the community may fetch and checkout your PR,
but it would usually be for testing out, or making minor suggestions.

Otherwise, if your commit is too old, or contains many redundant changes,
it would be a better practice to rebase and squash before asking for merge,
because currently, **squash-and-merge** is not enabled in SymPy project.

#### How to merge

First merge your local repository with the remote:

```
git checkout master
git pull
```

This results in:

```
A---B---C---D       master
         \
          a---b     fix-solve-bug
```

Then merge your `fix-solve-bug` branch from `fix-solve-bug`:

```
git checkout fix-solve-bug
git merge master
```

If the last command tells you that conflicts must be solved for a few indicated files.

If that's the case then the marks **>>>** and **\<\<\<** will appear at those files. Fix the
code with **>>>** and **\<\<\<** around it to what it should be.
You must manually remove useless pieces, and leave only new changes from your branch.

Then be sure that all tests pass:

```
./bin/test
./bin/doctest
```

and commit:

```
git commit
```

So the result will be like that (automatic merging `c`):

```
A---B---C-------D     master
         \       \
          a---b---M   fix-solve-bug
```

#### How to rebase

**If you have already made a pull request, please merge instead of rebasing.**

The final aim, that we want to obtain is:

```
A---B---C---D           master
             \
              a---b     fix-solve-bug
```

The way to do it is first of all to merge local repository with the remote `sympy/sympy`:

```
git checkout master
git pull
```

So we obtain:

```
A---B---C---D       master
         \
          a---b     fix-solve-bug
```

Then:

```
git checkout fix-solve-bug
git rebase master
```

Note that this last one will require you to fix some merge conflicts if there are changes
to the same file in `master` and `fix-solve-bug`. Open the file that it tells you is wrong,
fix the code with **>>>** and **\<\<\<** around it to what it should be.

Then be sure that all tests pass:

```
./bin/test
./bin/doctest
```

Then do:

```
git add sympy/matrices/your_conflict_file
git rebase --continue
```

(git rebase will also guide you in this).

#### Merging vs Rebasing

It is important to note that since rebase rewrites history, it is possible to
lose data, and it makes it harder for people reviewing your code, because they
can no longer just look at the "new commits"; they have to look at everything
again, because all the commits are effectively new.

There are several advantages to merging instead of rebasing. Rebasing
reapplies each commit iteratively over master, and if the state of the files
changed by that commit is different from when it was originally made, the
commit will change. This means what you can end up getting commits that are
broken, or commits that do not do what they say they do (because the changes
have been "rebased out"). This can lead to confusion if someone in the future
tries to test something by checking out commits from the history. Finally,
merge conflict resolutions can be more difficult with rebasing, because you
have to resolve the conflicts for each individual commit. With merging, you
only have to resolve the conflicts between the branches, not the commits.  It
is quite common for a merge to not have any conflicts but for a rebase to have
several, because the conflicts are "already resolved" by later commits.

Merging keeps everything intact. The commits you make are exactly the same,
down to the SHA1 hash, which means that if you checkout a commit from a merged
branch, it is exactly the same as checking it out from a non-merged branch.
What it does instead is create a single commit, the merge commit, that makes
it so that the history is both master and your branch.  This commit contains
all merge conflict resolution information, which is another advantage over
rebasing (all merge conflict resolutions when rebasing are "sifted" into the
commits that caused them, making them invisible).

#### Changing of commit messages

The only time when it is recommended to rebase instead of merge is when you
need to edit your commit messages, or remove unnecessary commits.

Note, it is much better to get your commit messages right the first time. See
the section on [writing good commit
messages](development-workflow-commit-messages) above.

Consider these commit messages:

```
$ git log --oneline
7bbbc06 More bug fixes
4d6137b Some additional corrections
925d88fx Fix a bug in solve()
```

Then run the *rebase* command in interactive mode:

```
git rebase --interactive 925d88fx
```

Or you can use other ways to point to commits, e.g. `git rebase --interactive HEAD^^`
or `git rebase --interactive HEAD~2`.

A new editor window will appear (note that order is reversed with respect to the `git log` command):

```
pick 4d6137b Some additional corrections
pick 7bbbc06 More bug fixes

# Rebase 925d88f..7bbbc06 onto 925d88f
#
# Commands:
#  p, pick = use commit
#  r, reword = use commit, but edit the commit message
#  e, edit = use commit, but stop for amending
#  s, squash = use commit, but meld into previous commit
#  f, fixup = like "squash", but discard this commit's log message
```

To edit a commit message, change *pick* to *reword* (or on old versions of
git, to *edit*) for those that you want to edit and save that file.

To squash two commits together, change *pick* to *squash*. To remove a commit,
just delete the line with the commit.

To edit a commit, change *pick* to *edit*.

After that, git will drop you back into your editor for every commit you want to reword,
and into the shell for every commit you wanted to edit:

```
git commit --amend -m "your new message"
git rebase --continue
```

For commits that you want to edit, it will stop. You can then do:

```
git reset --mixed HEAD^
```

This will "uncommit" all the changes from the commit. You can then recommit
them however you want. When you are done, remember to do:

```
git rebase --continue
```

Most of this sequence will be explained to you by the output of the various commands of git.
Continue until it says:

```
Successfully rebased and updated refs/heads/master.
```

If at any point you want to abort the rebase, do:

```
git rebase --abort
```

**Warning**: this will run `git reset --hard`, deleting any uncommitted
changes you have. If you want to save your uncommitted changes, run `git
stash` first, and then run `git stash pop` when you are done.


### Add your name and email address to the .mailmap file.

Every author's name and email address is stored in the AUTHORS file but this file should not be edited directly. The AUTHORS file is updated automatically when a new version of SymPy is released based on the name and email addresses that are recorded in the commits. Every commit made with git stores the name and email address that git is configured with (see "Configure git settings" above). The .mailmap file is used to associate the name/email recorded in the commits with an author name and email address that will be listed in the AUTHORS file.

The first time you make a pull request you will need to add your name and email address to the .mailmap file by adding a line like

```
Joe Bloggs <joe@bloggs.com>
```

This name and email should exactly match the name and email that you have configured with git before making the commits (see "Configure git settings" above). The `bin/mailmap_check.py` script can check that this has been done correctly. If you have made a commit but not yet added yourself to the .mailmap file then you will see this:

```bash
$ python bin/mailmap_check.py
This author is not included in the .mailmap file:
Joe Bloggs <joe@bloggs.com>

The .mailmap file needs to be updated because there are commits with
unrecognised author/email metadata.


For instructions on updating the .mailmap file see:
https://github.com/sympy/sympy/wiki/Development-workflow

The following authors will be added to the AUTHORS file at the
time of the next SymPy release.
```

This means that you should add your name and email address to the .mailmap file. If you add this at the end of the file then `git diff` will show:

```bash
$ git diff
diff --git a/.mailmap b/.mailmap
index 3af6dc1..7fa63b1 100644
--- a/.mailmap
+++ b/.mailmap
@@ -1307,3 +1307,4 @@ zsc347 <zsc347@gmail.com>
 Øyvind Jensen <jensen.oyvind@gmail.com>
 Łukasz Pankowski <lukpank@o2.pl>
 彭于斌 <1931127624@qq.com>
+Joe Bloggs <joe@bloggs.com>
```

Now you can rerun the `bin/mailmap_check.py` script and you should see:

```bash
$ python bin/mailmap_check.py
The mailmap file was reordered

For instructions on updating the .mailmap file see:
https://github.com/sympy/sympy/wiki/Development-workflow

The following authors will be added to the AUTHORS file at the
time of the next SymPy release.

Joe Bloggs <joe@bloggs.com>
```

The first line their says that the .mailmap file was "reordered". This is because the file should be in alphabetical order. The script will have moved your name into the correct position so now you can see the change as:

```bash
$ git diff
diff --git a/.mailmap b/.mailmap
index 3af6dc1..7598d94 100644
--- a/.mailmap
+++ b/.mailmap
@@ -562,6 +562,7 @@ Joannah Nanjekye <joannah.nanjekye@ibm.com> Joannah Nanjekye <jnanjekye@python.o
 Joannah Nanjekye <joannah.nanjekye@ibm.com> nanjekyejoannah <joannah.nanjekye@ibm.com>
 Joaquim Monserrat <qmonserrat@mailoo.org>
 Jochen Voss <voss@seehuhn.de>
+Joe Bloggs <joe@bloggs.com>
 Jogi Miglani <jmig5776@gmail.com> jmig5776 <jmig5776@gmail.com>
 Johan Blåbäck <johan_bluecreek@riseup.net> <johan.blaback@cea.fr>
 Johan Guzman <jguzm022@ucr.edu>
```

Now if you rerun the script you will see:

```bash
$ python bin/mailmap_check.py
No changes needed in .mailmap

The following authors will be added to the AUTHORS file at the
time of the next SymPy release.

Joe Bloggs <joe@bloggs.com>
```

The key information here is "No changes needed in .mailmap" which means that you have correctly updated the .mailmap file. You should now add and commit these changes as well:

```bash
git add .mailmap
git commit -m 'author: add Joe Bloggs to .mailmap'
```

Sometimes a commit will be made with an incorrect name or email address or an author will make multiple commits with different names and email addresses. In this case a line should be added to the .mailmap file where the first name and email address is what should be recorded in the AUTHORS file and the others are the name and email address that was incorrectly used in the other commits. For example if the commit was recorded with the name `joeb` and the email address `wrong@email.com` but the AUTHORS file should show `Joe Bloggs` as above then there should be a line in the .mailmap file like:

```
Joe Bloggs <joe@bloggs.com> joeb <wrong@email.com>
```

A common reason that this can happen is if making commits with the GitHub web UI which always recorded the name as github username and the email as something like `1785690389+joeb@users.noreply.github.com`. In this case a line will need to be added to .mailmap like:

```
Joe Bloggs <joe@bloggs.com> joeb <1785690389+joeb@users.noreply.github.com>
```

Multiple lines like this can be added to the .mailmap file. They should record all of the different name and email address combinations that have been used by an author and map all of them to a single author name that will show in the AUTHORS file.

If your pull request is merged and you have not previously been added to the AUTHORS file then your name will be added at the time of the next release of SymPy.

(development-workflow-reviewing)=
## Reviewing Pull Requests

Coding's only half the battle in software development: our code also has to be
thoroughly reviewed before release. Reviewers thus are an integral part of the
development process. Note that you do *not* have to have any special pull
or other privileges to review patches: anyone with Python on his/her computer
can review.

Every pull request to SymPy must be reviewed by someone else and signed-off on
before it can be merged. This applies equally to new contributors and to
existing maintainers of the project.

Feel free to view any open [pull
request](https://github.com/sympy/sympy/pulls) and discuss it. Reviews are the
place to read the changes and check them for correctness, and to ensure they
follow the [requirements for
inclusion](development-workflow-requirements-for-inclusion), such as including
tests.

Reviewing serves several purposes.

- To discuss technical implementation details and the best way to implement a
  fix or improvement.
- To ensure the changes are correct, both in terms of mathematics and the
  code logic.
- To ensure the changes are have a code quality that is up to the standards of
  the project (see [](development-workflow-requirements-for-inclusion)).
- To help the contributor with issues they may be having, either with a
  specific detail of the pull request, or with issues with git.

Some tips for reviewing:

- If you are requesting a change that is simple, use the "suggested change"
  feature in GitHub. This will allow the PR author to commit the change
  directly from the GitHub web interface. This is useful for things like typo
  and spelling fixes.

- Do not get too caught up in stylistic issues. It's OK if some code is not
  written in exactly the way you would have written it. At the same time, we
  do want SymPy's code to be high quality, so it is OK to request improvements
  if a piece of code is too low quality.

- Do not request changes for things that are unrelated to the main point of
  the pull request. Such things can be done in separate pull requests. Doing
  this tends to expand the scope of a pull request too much, and makes it
  harder to get merged, and leads to an endless review cycle where there are
  constantly new changes being requested. This also applies to new features
  that are not necessary for a pull request to be functional.

  [This blog
  post](https://matthewrocklin.com/blog/work/2022/03/26/fast-review) by
  Matthew Rocklin goes over why it's important to avoid scope creep when
  reviewing pull requests.

- If you have push access, you can use the feature that lets you push directly
  to someone's branch. This is useful if there is only a small change or two
  that needs to be made before the PR can be merged. This lets you push the
  change and merge it, without having to wait for the author. This is also
  useful for older pull requests which may be abandoned by the authors.

- Try to not let pull requests sit for too long without review. It is
  discouraging to contributor to make a pull request and have it sit
  unreviewed. However, as a contributor, you should also understand that
  SymPy's reviewers do so on a volunteer basis.


### Manual testing

If you prefer to test code manually, you will first have to set up your
environment as described in [](#create-your-environment). Then, you need to
obtain the patched files. If you're reviewing a pull request, you should get
the requested branch into your sympy folder. Go into your folder and execute
(\<username> being the username of the pull requester and \<branchname> being
the git branch of the pull request):

```
git remote add <username> git://github.com/<username>/sympy.git
git fetch <username>
git checkout -b <branchname> <username>/<branchname>
```

After obtaining the pull request or patch, go to your sympy root directory and
execute:

```
./bin/test
./bin/doctest
```

If there are any problems, notify the author in the pull request by commenting.

(development-workflow-requirements-for-inclusion)=
### Requirements for inclusion

A pull request or patch must meet the following requirements during review
before being considered as ready for release.

- All tests must pass.
  - Rationale: We need to make sure we're not releasing buggy code.
  - If new features are being implemented and/or bug fixes are added,
    tests should be added for them as well.
- The reviews (at least 1) must all be positive.
  - Rationale: We'd like everyone to agree on the merits of the patch.
  - If there are conflicting opinions, the reviewers should reach a consensus.
- The patch must have been posted for at least 24 hours.
  - Rationale: This gives a chance for everyone to look at the patch.

[doctests]: https://docs.python.org/3/library/doctest.html
[easy to fix issues]: https://github.com/sympy/sympy/issues?q=is%3Aopen+is%3Aissue+label%3A%22Easy+to+Fix%22
[GitHub]: https://github.com/
[Gitter sympy/sympy]: https://gitter.im/sympy/sympy
[issue tracker]: https://github.com/sympy/sympy/issues
[license]: https://github.com/sympy/sympy/blob/master/LICENSE
[mailing list]: https://groups.google.com/group/sympy
[new bsd license]: http://en.wikipedia.org/wiki/BSD_licenses#3-clause_license_.28.22New_BSD_License.22_or_.22Modified_BSD_License.22.29
[sympy]: https://sympy.org/
[sympy/sympy]: https://github.com/sympy/sympy
[sympy@googlegroups.com]: https://groups.google.com/group/sympy
