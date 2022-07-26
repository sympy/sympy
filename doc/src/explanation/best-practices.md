# Best Practices

This page outlines some of the best practices for users of SymPy. The best
practices here will help avoid some common bugs and pitfalls that can occur
when using SymPy.

## An Important Aside: Interactive vs. Programmatic Usage

There are two primary modes of usage for SymPy, and Python in general,
interactive and programmatic.

Interactive usage refers to live usage in an
interactive Python session such as a Python or [IPython](https://ipython.org/)
terminal session or [Jupyter Notebook](https://jupyter.org/). Interactive
sessions focus on live experimentation and fast iteration. It typically
involves reaching a solution in small iterative steps, viewing outputs often,
and going back and redefining parts of a calculation. Interactive sessions are
ephemeral in nature. The session history is either not saved at all, or even
if it is saved, it is not written or edited with any intention of rerunning
the code at an unspecified point in the future.

To contrast programmatic usage is any usage where the code saved with the
intention of being run again in the future. Typical interactive usage of SymPy
involves code that is saved to a `.py` file, or to a Jupyter notebook file
that is intended to be rerun. In programmatic usage, iterative steps are
consolidated to make code more readable, and abandoned ideas are deleted
entirely.

The distinction between interactive and programmatic usage matters because
many of the best practices outlined here may be ignored in interactive usage.
This is true of many general coding best practices, not just in SymPy. For
example, it is a best practice in Python to avoid wildcard imports, like `from
sympy import *`. These wildcard imports save typing, but they make it harder
for someone reading the code to understand where things come from, and they
can lead to bugs if wildcard imports from multiple libraries are used. But for
interactive usage, a wildcard import is acceptable, because it does save on
typing, and the ambiguity about what is actually imported can be resolved by
examining the namespace live in the interactive session.

SymPy has similar "shortcuts", which are perfectly acceptable in interactive
usage, but should be avoided for programmatic usage. Such instances will be
notated where appropriate. Other best practices listed here, which are not
notated as such, should be followed in all cases, as they can prevent pitfalls
ans bugs even in interactive usage.

## Basic Usage

### Defining Symbols

### Avoid String Inputs

### Exact Rational Numbers vs. Floats

## Custom SymPy Objects

### Args Invariants

### Avoiding Type Confusion

#### Don't Denest Collections
