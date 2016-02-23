# Prepare the VM

First execute:

    vagrant up
    fab vagrant prepare

which will prepare the VM (install packages, cache sympy repository, etc.).

You only need to execute this once. It will take a while if you have never run
it before, because it has to download a lot of stuff.

# Release

First, make sure that you have done the following things

- Create a release branch. Usually this branch is the same name as the release
  (e.g., "0.7.3"), although no naming convention is enforced on it.

- Change the version in the release branch in sympy/release.py. If you want to
  do a release candidate, change it to a [PEP
  440](https://www.python.org/dev/peps/pep-0440) compliant version like
  0.7.3.rc1.

- Change the version in master. This way, any additional changes made in master
  will be shown as coming from the right place. The master release should be
  e.g. `0.7.4.dev`, see [PEP 440](https://www.python.org/dev/peps/pep-0440) for
  rules about development version numbers. Note that this version number should
  the next projected version plus the `.dev`.

- Push the release branch up to origin, and make a pull request for it against
  master.

It is important to create a new branch because that lets master continue
as normal. The fab script will automatically checkout the release branch from
origin, which is why you need to push it (it determines what the release
branch by just looking at what branch you have checked out locally, so make
sure you are on the release branch when you release). It is important to
change the version number because it uses that in naming the tarballs it
creates.

If you want to test the release process without pushing a branch to the
official repo, you can push a branch to your fork and use `fab vagrant
release:fork='username'`, where `username` is your GitHub username.  Note that
once you do the actual release, you should do it in a branch in the official
GitHub repo. **NOTE**: If your fork does not have all the tags of the
official repo, then the code that finds the previous version will not work
correctly.  Hence, you may see things like more authors in the authors list
than you should.  To remedy this, be sure to do `git fetch origin --tags` and
`git push github --tags`.

Once you have done these things, execute:

    fab vagrant release

this create release tarballs and put them all into a new "release" directory
of the current directory.

# Testing things

The full test suite is not run by fabric, because we leave that to
Travis. However, there are things that need to be tested specific to the
release. Most of these things are done automatically by the release command
(like testing that the tarball can be installed), but one thing must be tested
manually, because it has to be inspected by hand, namely, making sure that the
tarballs contain everything, and don't contain any junk files.

Run

    fab vagrant show_files:arg

to show the files in the tarball, where `arg` is `source` or `html`.  You'll
probably want to pipe the output of this into `less`, so that you can inspect
it.

You should also open the pdf and make sure that it has built correctly, and
open the html docs and make sure that they have built correctly.

# Tagging the release

Once you have made the final release files that you plan to upload, be sure
that everything is committed, and that the most recent git HEAD is indeed the
same one that was used to build the files (you can always do `fab vagrant
release` again if you are not sure). Then tag the release with the command

    git tag sympy-VERSION -a

where you should replace `VERSION` with the version (which should be `x.y.z`,
or `x.y.z.rcn` for the `n`th release candidate. It is very important to follow
the tag naming conventions.  The `-a` will cause it to prompt for a tag commit
message. Just write something like "SymPy VERSION release".

Then, push up the tag, with

    git push origin sympy-VERSION

Note, once a tag is pushed, that's it. It can't be changed. If you need to
change the tag, you must bump the release number.  So double check that
everything is right before pushing.

# Uploading

Before you release, you need to push the tag up, as described above.

Release candidates should be uploaded to GitHub only.

    fab vagrant GitHub_release

This will create the release on GitHub for the tag, and upload the files to
it.  Do not upload release candidates to PyPI, as `pip` and `easy_install`
will pick them up if you do.

This will prompt you for a username and password the first time you call it.
After that, it will prompt you to generate a token file.  If you don't save
the token to a file, you will need to pass it in as an argument. Releasing is
only supported via OAuth, so using a token is required.

You (obviously) need push access to create a GitHub release.

If you want to test this before doing it, use

    fab vagrant GitHub_release:draft=True

This will make the release not visible until you go to the web interface and
publish it.  You can also set the `user` and `repo` flags to test against a
different GitHub repo.

For final releases, you should upload to both GitHub and PyPI. The command

    fab vagrant upload

will do both of these.  You will need admin access to the SymPy PyPI project.

Note that if either of these commands fails for some reason, you will very
likely need to go into the web interface and clean some things up before you
can upload again.

# Updating websites

You should now update the websites. Only do this for final releases. The command

    fab vagrant update_websites

will update docs.sympy.org and sympy.org.  You will need to have local clones
of these repos, and push access to them (obviously).  **Note, this command
will commit and push the changes automatically.**

The other website that needs to be updated is SymPy Live. You should make this
as a pull request to the Live repo.

# Other

You can run all the SymPy tests by running:

    fab vagrant test_sympy

To get the md5 sums of all the files, use

    fab md5

To list the files in the tarball use

    fab vagrant show_files:arg

where `arg` is `source` or `html`, for the Python sources and the html docs,
respectively. Note that the source code is already checked automatically
against the files in git and a whitelist.

You can obtain all available commands by:

    fab -l

# Restarting from scratch

Run

    vagrant destroy

You can also delete the releases that it has built

    rm -rf release
