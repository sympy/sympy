.. _devsetup:

================================
Development Environment Setup
================================

The first step to contributing to the code base is creating your development environment.

Git Setup
-----------

SymPy is available on `GitHub <https://github.com/sympy/sympy>`_ and uses
`Git <https://git-scm.com>`_ for source control. The workflow is such that
code is pulled and pushed to and from the main repository. Install the respective version
of Git for your operating system to start development.

.. note::
   Refer to the installation instructions in
   the `Git installation instructions <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_.
   Learn about the basic git commands in this `Git Handbook <https://guides.github.com/introduction/git-handbook/>`_
   or any other sources on the internet.

Get the SymPy Code
-------------------

It is recommended practice to create a fork of the SymPy project for your development purposes. Create your own fork of the SymPy project (if you have not yet). Go to the SymPy GitHub repository:

.. code-block:: bash

   https://github.com/sympy/sympy

You will now have a fork at <https://github.com/<your-user-name>/sympy>.

Then, nn your machine browse to where you would like to store SymPy, and clone (download) the latest code from SymPy's original repository (about 77 MiB):

.. code-block:: bash

   $ git clone https://github.com/<your-user-name>/sympy

You must `configure the remote repositories <https://git-scm.com/book/en/v2/Git-Basics-Working-with-Remotes>`_ for collaboration with the upstream project:

.. code-block:: bash

   $ cd sympy
   $ git remote add upstream https://github.com/sympy/sympy

After the configuration, your setup should be similar to this:

.. code-block:: bash

   $ git remote -v
   origin   https://github.com/<your-user-name>/sympy (fetch)
   origin   https://github.com/<your-user-name>/sympy (push)
   upstream https://github.com/sympy/sympy (fetch)
   upstream https://github.com/sympy/sympy (push)

For further development, it is recommended
to create a development branch.

.. code-block:: bash

    $ git checkout -b dev-branch

The new branch can be of any name.

Virtual Environment Setup
---------------------------

You may want to take advantage of using virtual environments to isolate your development version of SymPy from any system wide installed versions, e.g. from ``apt-get install python-sympy``.

We recommend using ``conda`` to create a virtual environment:

.. code-block:: bash

    $ conda create -n sympy-dev python=3 mpmath flake8

You now have a environment that you can use for testing your development copy of SymPy. For example, clone your SymPy fork from Github:

.. code-block:: bash

    $ git clone git@github.com:<your-github-username>/sympy.git
    $ cd sympy

Now activate the environment:

.. code-block:: bash

    $ conda activate sympy-dev


Run the Tests
--------------

There are several ways of running SymPy tests but the easiest is to use the ``bin/test`` script, consult `the wiki details on running tests <https://github.com/sympy/sympy/wiki/Running-tests>`_.

The script takes a number of options and arguments and then passes them to ``sympy.test(*paths, **kwargs)``. Run ``bin/test --help`` for all supported arguments.

Run all tests by using the command:

.. code-block:: bash

    $ bin/test

To run tests for a specific file, use:

.. code-block:: bash

    $ bin/test test_basic

Where ``test_basic`` is from file ``sympy/core/basic.py``.

To run tests for modules, use:

.. code-block:: bash

   $  bin/test /core /utilities

This will run tests for the ``core`` and ``utilities`` modules.

Similarly, run quality tests with:

.. code-block:: bash

    $ bin/test code_quality
