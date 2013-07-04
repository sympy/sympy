# Prepare the VM

First execute:

    vagrant up
    fab vagrant prepare

which will checkout SymPy and prepare the VM.

# Release

Execute:

    fab vagrant release

# Other

You can test SymPy by:

    fab vagrant sympy_test

You can obtain all available commands by:

    fab -l
