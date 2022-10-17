SymPy is created and maintained by a large group of contributors and we'd love you to be one of them too! Getting started in a well-oiled large and complex machine like SymPy can be daunting for new contributors. This page exists to give tips for new contributors to get started.

Get familiar using the software
================================

We suggest going through the [SymPy Tutorial](http://docs.sympy.org/latest/tutorial/) to get acquainted with using the software before you begin contributing. This will help you to familiarize yourself with the uses of SymPy.

The tutorial is also available on video:

- [Symbolic Computation with Python using SymPy | SciPy 2016](https://www.youtube.com/watch?v=AqnpuGbM6-Q)
- SymPy Tutorial SciPy 2014 [Part 1](https://www.youtube.com/watch?v=Lgp442bibDM) [Part 2](https://www.youtube.com/watch?v=_PTe10whFKo) [Part 3](https://www.youtube.com/watch?v=qleGSnrnxgc)

Read the paper
==============

We authored a journal paper in 2017 that provides a high-level look at SymPy and its capabilities. You can read it here:

https://peerj.com/articles/cs-103/

Peruse the documentation
========================

Besides the tutorial, there is a lot more information in the [documentation](http://docs.sympy.org/latest/index.html). It's probably a good idea to at least browse through the different topics to get an idea of what else is available. For example, the [architecture](https://web.archive.org/web/20200403130424/https://docs.sympy.org/latest/guide.html) section gives an overview of how the SymPy code is organized and is useful for a contributor to know about.

Review the Code of Conduct
==========================

Participants in the SymPy community are expected to abide by our [Code of Conduct](https://github.com/sympy/sympy/blob/master/CODE_OF_CONDUCT.md). Please review this before getting started.

Join our mailing list
=====================

The [SymPy email mailing list](https://groups.google.com/forum/#!forum/sympy) is one place where discussions about SymPy happen. You can ask questions about how to use SymPy, discuss feature requests, discuss software bugs, or share how you are using SymPy. Request to join the list on the Google Groups page. Please read http://shakthimaan.com/downloads/book/chapter1.pdf before posting to get familiar with mailing list etiquette.

If posting for the first time then it is good to send an informative first email so that we can get to know you. Here are some example pieces of information that would be useful for you to include in an introductory email:

- level of familiarity with python (years of programming, previous projects)
- mathematical education level (high school / ... / PhD?)
- particular expertise (physics? biology? ...subtopics)
- particular algorithmic interests
- level of familiarity with symbolic math systems "computer algebra". Perhaps
even course work in (say) algebra or books read.
- your familiarity with SymPy (e.g. how have you used SymPy)?
- other possibly relevant information -- geographical location? native language?

Join the Gitter chat room
=========================

We have a live chat room where we discuss things about SymPy at Gitter: https://gitter.im/sympy/sympy. Here you may find other SymPy developers and you can ask questions in real-time.

Peruse the Wiki
===============

Another good source of information is the [SymPy Wiki](https://github.com/sympy/sympy/wiki). It's a good idea to browse through the topics there. And feel free to create or edit a page. For example, you can create a page to plan out code modifications you have in mind.

Setup your development environment
==================================

We use the [Git](http://git-scm.com) version control system to track the software [changes over time](https://github.com/sympy/sympy/commits/master) and to effectively manage [contributions from many different authors](https://github.com/sympy/sympy/network). We also utilize Github, a web interface to Git, extensively and use it for communication, issue tracking, merging patches (i.e. Pull Requests), etc.

You will need to read through the [[Development Workflow]] page and follow the instructions to set up your development environment before making a contribution. It is easiest for us if you make use of Github's "Pull Request" system for sending in patches and code contributions.

Identify something to work on
=============================

There are lots of ways to contribute to SymPy. Most contributions center around fixing software bugs and adding new features. But there are other things we need help with too, like maintaining our websites, writing documentation, preparing tutorials, answering people's questions on the mailing lists, chat room, and issue tracker. Here are some following ways to get started with a contribution:

SymPy Codebase
--------------

The best way to start with the main codebase is to fix some existing bugs. Peruse the ["Easy to fix" issues](https://github.com/sympy/sympy/issues?q=is%3Aopen+is%3Aissue+label%3A%22Easy+to+Fix%22) in the issue tracker and see if one interests you. If you'd like to try to fix it, then create a message in the issue saying that you'd like to work on it. If it isn't clear how to fix it, ask for suggestions on how to do it in the issue itself, on the mailing list, or on Gitter.

SymPy's code is organized into Python packages and modules. The core code is in the `/sympy/core` directory and other packages in the `/sympy` directory have more specific code, for example `/sympy/printing` handles how SymPy objects are printed to the terminal, in IPython notebooks, or in our web applications.

Project Ideas Page
------------------

If you are looking for a somewhat larger project to implement, check out the [Project General Ideas](https://github.com/sympy/sympy/wiki/Project-General-Ideas) page. This page is a collection of projects that contributors have come up with but have not yet had the time or opportunity to implement themselves. Another good place to look for ideas are the Google Summer of Code ideas pages: [[GSoC 2018 Ideas]], [[GSoC 2017 Ideas]], [[GSoC 2016 Ideas]], ...

Documentation
-------------

If you'd like to improve the documentation you can edit in one of two places:

1. The documentation source files: https://github.com/sympy/sympy/tree/master/doc/src
2. The docstrings* of the functions in the source code: https://github.com/sympy/sympy/tree/master/sympy

Both of these end up displayed on the documentation website:

1. Example rendered prose: https://docs.sympy.org/dev/guides/index.html and it's corresponding source: https://github.com/sympy/sympy/tree/master/doc/src/guides
2. Example rendered docstrings: http://docs.sympy.org/latest/modules/core.html and the corresponding source: https://github.com/sympy/sympy/blob/master/sympy/core/sympify.py#L55

\* Every function and class in SymPy has a string below call signature explaining the use of the object. This is what is displayed in Python when you type `help(function_name)`.

While contributing to or improving upon our documentation, please follow the [SymPy Documentation Style Guide](https://docs.sympy.org/dev/documentation-style-guide.html).

SymPy Websites
--------------

We have three websites: [sympy.org](http://sympy.org/en/index.html), [live.sympy.org](http://live.sympy.org/), [sympygamma.com](http://www.sympygamma.com/). The first is our primary web presence and is a "static" website. The last two are interactive web applications where you can use SymPy without installing it and they run on servers provided by the Google App Engine.

1. Edit the files for the main SymPy website: https://github.com/sympy/sympy.github.com and check out the issue tracker for the website: https://github.com/sympy/sympy.github.com/issues.
2. Both SymPy Live and SymPy Gamma have their own repositories: https://github.com/sympy/sympy-live https://github.com/sympy/sympy_gamma and issue trackers. If you'd like to contribute to those websites, start there.

Review pull requests
--------------------

Every repository has a "Pull Request" section where people send the code they'd like to contribute to SymPy. For example, https://github.com/sympy/sympy/pulls. You can view the code submission and check whether it does what it is intended to do. Every pull request has to be reviewed and approved by a reviewer before we merge it into the main code base. The blogpost ["Mindful Communications in Code Reviews"](https://kickstarter.engineering/a-guide-to-mindful-communication-in-code-reviews-48aab5282e5e) is helpful to read before reviewing for SymPy. 