# Prepare the VM

First execute:

    vagrant up
    fab vagrant prepare:0.7.3

which will checkout SymPy into the branch 0.7.3 (you can use any branch name)
and prepare the VM.

# Release

Execute:

    fab vagrant release

# Other

You can test SymPy by:

    fab vagrant sympy_test

You can obtain all available commands by:

    fab -l
