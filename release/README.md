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
(e.g., "0.7.3").

- Change the version in the release branch in sympy/__init__.py.  If you want
  to do a release candidate, change it to something like 0.7.3.rc1.

- Change the version in master.  This way, any additional changes made in
  master will be shown as coming from the right place. The master release
  should be like "0.7.3-git".

- Push the release branch up to origin, and make a pull request for it against
  master.

It is important to create a new branch because that lets master continue
as normal. The fab script will automatically checkout the release branch from
origin, which is why you need to push it (it determines what the release
branch by just looking at what branch you have checked out locally, so make
sure you are on the release branch when you release). It is important to
change the version number because it uses that in naming the tarballs it
creates.

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

to show the files in the tarball, where `arg` is one of `2`, `3`, or `html`.
You'll probably want to pipe the output of this into less, so that you can
inspect it.

You should also open the pdf and make sure that it has built correctly, and
open the html docs and make sure that they have built correctly.

# Uploading

TODO

# Other

You can run all the SymPy tests by running:

    fab vagrant test_sympy

To get the md5 sums of all the files, use

    fab md5

You can obtain all available commands by:

    fab -l

# Restarting from scratch

Run

    vagrant destroy

You can also delete the releases that it has built

    rm -rf release
