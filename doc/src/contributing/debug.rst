
Debugging
==========

To start sympy in debug mode set the SYMPY_DEBUG variable. For instance in a unix-like system you would do

    $ SYMPY_DEBUG=True bin/isympy

or in Windows

    > set SYMPY_DEBUG=True
    > python bin/isympy

Now just use for example the ``limit()`` function. You will get a nice printed tree, which is very useful for debugging.
