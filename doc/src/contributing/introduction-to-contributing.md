# Introduction to Contributing

SymPy is created and maintained by a large group of contributors and we'd love you to be one of them too! Getting started in a well-oiled large and complex machine like SymPy can be daunting for new contributors. This page exists to give tips for new contributors to get started.

## Get familiar using the software

We suggest going through the [SymPy Tutorial](intro-tutorial) to get acquainted with using the software before you begin contributing. This will help you to familiarize yourself with the uses of SymPy.

The tutorial is also available on video:

- [Symbolic Computation with Python using SymPy | SciPy 2016](https://www.youtube.com/watch?v=AqnpuGbM6-Q)
- SymPy Tutorial SciPy 2014 [Part 1](https://www.youtube.com/watch?v=Lgp442bibDM) [Part 2](https://www.youtube.com/watch?v=_PTe10whFKo) [Part 3](https://www.youtube.com/watch?v=qleGSnrnxgc)

## Read the paper

We authored a journal paper in 2017 that provides a high-level look at SymPy and its capabilities. You can read it here:

https://peerj.com/articles/cs-103/

## Peruse the documentation

Besides the tutorial, there is a lot more information in the [documentation](documentation). It's probably a good idea to at least browse through the different topics to get an idea of what else is available.

## Review the Code of Conduct

Participants in the SymPy community are expected to abide by our [Code of Conduct](https://github.com/sympy/sympy/blob/master/CODE_OF_CONDUCT.md). Please review this before getting started.

## Join our mailing list

The [SymPy email mailing list](https://groups.google.com/forum/#!forum/sympy)
is one place where discussions about SymPy happen. You can ask questions about
how to use SymPy, discuss feature requests, discuss software bugs, or share
how you are using SymPy. Request to join the list on the Google Groups page.
Note that to prevent spam, the first time you post your message will need to
be moderated before it is posted to the list. Please read
http://shakthimaan.com/downloads/book/chapter1.pdf before posting to get
familiar with mailing list etiquette.

## Setup your development environment

We use the [Git](https://git-scm.com) version control system to track the
software [changes over time](https://github.com/sympy/sympy/commits/master)
and to effectively manage [contributions from many different
authors](https://github.com/sympy/sympy/network). We also utilize GitHub, a
web interface to Git, extensively and use it for communication, issue
tracking, merging contributions (i.e., pull requests), etc.

If you are new to Git and GitHub, read through the [](devsetup) page first for
instructions on how to set up your development environment. If you are already
familiar with the basic GitHub workflow, the [](workflow-process) page
describes the aspects of the GitHub contribution workflow that are specific to
SymPy.

## Identify something to work on

There are lots of ways to contribute to SymPy. Most contributions center
around fixing software bugs and adding new features for things that are
interesting to them. But there are other things we need help with too, like
maintaining our websites, writing documentation, preparing tutorials,
answering people's questions on the mailing lists, chat room, StackOverflow,
and issue tracker, and reviewing pull requests. Here are some following ways
to get started with a contribution:

### SymPy Codebase

The best way to start with the main codebase is to fix some existing bugs. If
you are looking for a bug to fix, you can start by looking at the issues
labeled ["Easy to
fix"](https://github.com/sympy/sympy/issues?q=is%3Aopen+is%3Aissue+label%3A%22Easy+to+Fix%22)
in the issue tracker and see if one interests you. If it isn't clear how to
fix it, ask for suggestions on how to do it in the issue itself or on the
mailing list.

SymPy's code is organized into Python packages and modules. The core code is
in the `sympy/core` directory and other packages in the `sympy` directory have
more specific code, for example `sympy/printing` handles how SymPy objects are
printed to the terminal and Jupyter.

### Documentation

SymPy's documentation lives in two places:

1. The documentation source files: https://github.com/sympy/sympy/tree/master/doc/src
2. The docstrings* of the functions in the source code: https://github.com/sympy/sympy/tree/master/sympy

Both of these end up displayed here on the documentation website. You can
click the "[Source]" link next to any function documentation to go to its
corresponding docstring in the SymPy source code.

\* Every function and class in SymPy has a string below call signature explaining the use of the object. This is what is displayed in Python when you type `help(function_name)`.

While contributing to or improving upon our documentation, please follow the [SymPy Documentation Style Guide](documentation-style-guide).

### Review pull requests

Every contribution to SymPy goes through a pull request
https://github.com/sympy/sympy/pulls. We require that every pull request
undergo a review before it can be merged. If you have gained some familiarity
with a part of the SymPy codebase and the SymPy development processes, it can
be helpful to the community to review others' pull requests. You can view the
code submission and check whether it does what it is intended to do.
