===============
 Preliminaries
===============

This tutorial assumes that the reader already knows the basics of the Python programming
language.  If you do not, the `official Python
tutorial <http://docs.python.org/3/tutorial/index.html>`_ is excellent.

This tutorial assumes a decent mathematical background.  Most examples require
knowledge up to calculus level, and some require knowledge at a calculus
level.  Some of the advanced features require more than this. If you come
across a section that uses some mathematical function you are not familiar
with, you can probably skip over it, or replace it with a similar one that you
are more familiar with.  Or look up the function on Wikipedia and learn
something new.

You will need to install SymPy first.  See :ref:`installation`.

Alternately, you can just use the SymPy Live Sphinx extension to run the code
blocks in the browser.  For example, click on the green "Run code block in
SymPy Live" button below

    >>> from sympy import *
    >>> x = symbols('x')
    >>> a = Integral(cos(x)*exp(x), x)
    >>> Eq(a, a.doit())
    Integral(exp(x)*cos(x), x) == exp(x)*sin(x)/2 + exp(x)*cos(x)/2

The SymPy Live shell in the bottom corner will pop up and evaluate the code
block. You can also click an individual line to evaluate lines one at a time.

The SymPy Live shell is a fully interactive Python shell. You can type any
expression in the input box to evaluate it.  Feel free to use it throughout
the tutorial to experiment.

To show or hide the SymPy Live shell at any time, just click the green button
on the bottom right of the screen.

If you wish to modify an example before evaluating it, change the evaluation
mode to "copy".  This will cause clicking on an example to copy the example to
the SymPy Live shell, but not evaluate it, allowing you to change it before
execution.  You can also use the up/down arrow keys on your keyboard in the
input box to move through the shell history.

The SymPy Live shell is also available at http://live.sympy.org, with extra
features like a mobile phone enhanced version.
