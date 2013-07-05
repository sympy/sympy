# Prepare the VM

First execute:

    vagrant up
    fab vagrant prepare

which will prepare the VM (install packages, cache sympy repository, etc.).

# Release

Execute:

    fab vagrant release:0.7.3

this will checkout the branch that you specify (0.7.3 in this case), create
release tarballs and put them all into a new "release" directory of the current
directory.

# Other

You can test SymPy by:

    fab vagrant sympy_test

You can obtain all available commands by:

    fab -l
