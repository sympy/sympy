Release instructions
--------------------

Here we will describe the process of releasing version 1.12. At the start of
the release process the master branch records version `1.12.dev` and the
previous release was `1.11`. We will make a branch for the `1.12.*` release
series and the master branch will become `1.13.dev`.

1.  Make an issue to track the release. Link back to this issue from all PRs
    and issues related to the release. Use the OP to make a todo list. The
    release issue for 1.12 is here:
    https://github.com/sympy/sympy/issues/24601

2.  Announce the start of the release process on the mailing list.

3.  Clean up release blocker issues, especially regressions.

4.  Defer any issues that cannot be fixed for the current release.

5.  Update the AUTHORS file on the master branch *before* creating the release
    branch so that the update applies to both branches.

    1. Make sure you fully fetch the upstream master branch.
    2. Run `bin/mailmap_check.py`.
    3. Check for duplicates and fix the `.mailmap` file.
    4. Run `bin/mailmap_check.py --update-authors`.
    5. Open the PR and potentially use it to ask contributors to clarify what
       name they want to have recorded by `@`-tagging them in the PR.
    6. Check to see if any names in the AUTHORS file look like obvious
       duplicates and fix the .mailmap. Make sure to notify the authors and
       comment on the associated PR to make it clear that .mailmap has not been
       correctly edited.
    7. Most likely the AUTHORS file will not need to be updated on the release
       branch but if anyone does submit a PR to the release branch who hasn't
       submitted any previous PRs to SymPy then this will need to be done again
       on the release branch and ported to master as well.

    Only merge this PR just *before* creating the release branch (see below).

    Example from releasing SymPy 1.12:
    https://github.com/sympy/sympy/pull/24775

6.  Create a separate PR to update the version on the master branch from
    1.12.dev to 1.13.dev. The versions should be updated in:

    1. `sympy/release.py`
    2. `.github/workflows/release.yml`
    3. `doc/src/active-deprecations`

    Don't merge the PR before creating the release notes page for 1.12 because
    once the PR is merged any new PRs that are merged to master will have their
    release note go missing.

    Example from releasing SymPy 1.12:
    https://github.com/sympy/sympy/pull/24774

7.  Create the release notes page for 1.13. The release notes are automatically
    updated by the SymPy release notes bot after each PR is merged so it is
    important that the name and structure of the file is exactly correct. The
    release note page needs to exist before updating the version on the master
    branch.

    The release notes for 1.12 can be seen here:

    https://github.com/sympy/sympy/wiki/Release-Notes-for-1.12

    The title of the page should be `Release Notes for 1.13`.

    Use this template but change the versions of SymPy and check that the
    listed Python versions are up to date with what is tested in CI:

    ```markdown
    These are the release notes for SymPy 1.13. You can also find release notes for
    [[previous versions|https://github.com/sympy/sympy/wiki/Release-Notes]].

    SymPy 1.12 has not been released yet.

    This version of SymPy has been tested on Python 3.8, 3.9, 3.10, 3.11, and
    PyPy. See our [Python version support
    policy](https://github.com/sympy/sympy/wiki/Python-version-support-policy) for
    more information on when we plan to drop support for older Python versions.

    Install SymPy with

        pip install -U sympy

    or if you use Anaconda

        conda install sympy

    ## Highlights

    There are many changes in 1.12 (see below).

    ## Backwards compatibility breaks and deprecations


    ## Changes


    ## Authors
    ```
    Also add a link to the release notes page here:
    https://github.com/sympy/sympy/wiki/Release-Notes

8.  Create the release branch. These 3 steps need to happen immediately one
    after another in this order:

    1.  Merge the update to the AUTHORS file and check that no new authors have
        appeared since the PR was opened. There should be no output from
        `bin/mailmap_check.py` on the master branch after merging this.
    2.  Create the release branch and push to the main repo:
    ```console
    $ git fetch upstream
    $ git checkout upstream/master
    $ git checkout -b 1.12
    $ bin/mailmap_check.py  # Last check that AUTHORS and .mailmap are good.
    $ git push upstream 1.12
    ```
    3.  Merge the version update PR to the master branch. This needs to happen
        after creating the release branch so that no PRs are merged in between
        the release branch being created and the master branch version being
        updated (the bot might post the release note to the wrong page).

    If during the above you find that there are new .mailmap entries then
    create a new PR to update AUTHORS, delete the local 1.12 branch (without
    pishing to upstream) and go back to the first step). During this time any
    PRs that are merged will have their release notes go the 1.13 release notes
    page and you will need to move them over manually.

9.  Create a PR to update the versions on the release branch. This needs to
    update the same places.

    1. `sympy/release.py`
    2. `.github/workflows/release.yml`
    3. `.github/workflows/runtests.yml`
    4. `asv.conf.actions.json`

    In `runtests.yml` and `asv.conv.actions.json` the master branch should be
    replaced with the release branch (`1.12`).

10. Make a new PR to update the SymPy versions used for comparions in the
    benchmarks job in CI and also maybe to update Python versions used for
    running the tests. These files need to be updated:

    1. `.github/workflows/runtests.yml`
    2. `asv.conf.actions.json`

    We do this here rather than in the previous version bump PR because we
    needed to create the release branch before merging the changes for `asv`
    which will refer to the 1.12 release branch for benchmark comparisons.

    Example from releasing SymPy 1.12:
    https://github.com/sympy/sympy/pull/24776

11. Explain release process on mailing list and in release issue

    From this point on any new PRs merged to master will not go into the
    release by default. To go into the release they need to be backported to
    the release 1.11 branch. This should mostly be only bugfixes with a
    priority for any bug that is a regression introduced since the last release
    or in the last release.

    Various issues and PRs will have been marked with the 1.11 milestone
    implying that they are release blockers. Many of them will not really be
    release blockers though (i.e. they are not really things that can't just
    wait until the next release). Contributors should hurry up and finish any
    PRs that they want in the release but if they don't then the release
    manager can defer anything that is not a real blocker or just remove the
    milestone.
