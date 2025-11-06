=============================
Performance Benchmarking (ASV)
=============================

SymPy uses `Airspeed Velocity (ASV)`_ to measure and track performance
across commits. ASV allows contributors to identify regressions, compare
performance between branches, and validate the impact of their changes.

This guide explains how to install ASV, run benchmarks locally, interpret
results, and compare performance between commits.

.. contents::
   :local:
   :depth: 2

-----------------------
1. Installing and Setup
-----------------------

ASV is not installed automatically with SymPy. To install it:

.. code-block:: bash

    pip install asv

The benchmark suite is located in the top-level ``benchmarks/`` directory.

Before running ASV, ensure that SymPy is built in-place:

.. code-block:: bash

    python setup.py build_ext --inplace

This step is required so that ASV uses the local development version of SymPy.

--------------------------------
2. Running Benchmarks (Basic Use)
--------------------------------

To run *all* benchmarks on the current commit:

.. code-block:: bash

    asv run

This will run the suite according to the configuration in ``asv.conf.json``.

To run benchmarks more quickly (fewer iterations, fewer warmups):

.. code-block:: bash

    asv run --quick

To run benchmarks in a specific subdirectory:

.. code-block:: bash

    asv run benchmarks/solvers

Or a specific file:

.. code-block:: bash

    asv run benchmarks/solvers/bench_solve.py

To run only benchmarks whose names match a pattern:

.. code-block:: bash

    asv run --bench <pattern>

Example:

.. code-block:: bash

    asv run --bench integrate

------------------------------------
3. Comparing Results Across Commits
------------------------------------

ASV allows direct comparison between two revisions:

.. code-block:: bash

    asv compare master HEAD

The typical output looks like:

::

    benchmark                                      ratio    std
    ------------------------------------------------------------
    solve.SymbolicSolve.time_solve_quadratic       1.20     0.03

Interpretation:

* **ratio > 1.0** → the new commit is *slower* (possible regression)
* **ratio < 1.0** → the new commit is *faster*
* **std** shows timing variability (lower is more stable)

You can also compare arbitrary commits:

.. code-block:: bash

    asv compare <old> <new>

-----------------------------------
4. Interpreting ASV Benchmark Output
-----------------------------------

Each benchmark typically reports:

* **time** – execution time in seconds
* **ratio** – comparison factor between two commits
* **std** – standard deviation across runs
* **samples** – number of timing samples collected

General tips:

* Treat **ratio ≥ 1.2** as a likely regression.
* Large **std** often means system noise; rerun with ``asv run`` (without ``--quick``).
* For unstable benchmarks, running more samples may help.

-----------------------------------
5. Viewing Results with ASV Preview
-----------------------------------

To see interactive visualizations:

.. code-block:: bash

    asv preview

This opens a local web interface that includes:

* commit-by-commit performance graphs
* historical timing trends
* regression highlights
* benchmark descriptions
* environment configuration

This is recommended before submitting a pull request.

------------------------------------------------------
6. Detecting and Reporting Regressions in Pull Requests
------------------------------------------------------

Before creating a PR, contributors should check that their changes do not introduce performance regressions.

Recommended workflow:

1. Build SymPy in-place:

   .. code-block:: bash

       python setup.py build_ext --inplace

2. Run relevant benchmarks:

   .. code-block:: bash

       asv run benchmarks/<module>

3. Compare results with ``master``:

   .. code-block:: bash

       asv compare master HEAD

If a regression is found:

* Investigate whether the slowdown is expected or necessary.
* If unexpected, try profiling or isolating the cause.
* Mention ASV results in the PR description, including benchmark names and ratios.

Maintainers may request additional ASV comparisons before merging.

----------------------
7. Useful References
----------------------

* ASV documentation: https://asv.readthedocs.io/en/stable/
* SymPy benchmark suite: ``benchmarks/`` directory
* Contributor Guide: :ref:`contributing`

.. _Airspeed Velocity (ASV): https://asv.readthedocs.io/en/stable/