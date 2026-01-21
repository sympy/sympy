=============================
Performance Benchmarking (ASV)
=============================

SymPy uses Airspeed Velocity (ASV) to measure and track performance
across commits. ASV allows contributors to identify regressions, compare
performance between branches, and validate the impact of their changes.

This guide explains how to install ASV, run benchmarks locally, interpret
results, and compare performance between commits.

.. contents::
   :local:

----------------------
Installing and Setup
----------------------

ASV is not installed automatically with SymPy. To install it:

.. code-block:: bash

    pip install asv

The benchmark suite is located in the top-level ``benchmarks/`` directory.

Before running ASV, ensure that SymPy is built in-place:

.. code-block:: bash

    python setup.py build_ext --inplace

This step is required so that ASV uses the local development version of SymPy.

------------------------------
Running Benchmarks (Basic Use)
------------------------------

To run all benchmarks on the current commit:

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

----------------------------------
Comparing Results Across Commits
----------------------------------

ASV allows direct comparison between two revisions:

.. code-block:: bash

    asv compare master HEAD

Example output:

::

    benchmark                                      ratio    std
    ------------------------------------------------------------
    solve.SymbolicSolve.time_solve_quadratic       1.20     0.03

Interpretation:

* **ratio > 1.0** – slower (possible regression)
* **ratio < 1.0** – faster
* **std** – timing variability (lower is better)

To compare arbitrary commits:

.. code-block:: bash

    asv compare <old> <new>

----------------------------------
Interpreting ASV Benchmark Output
----------------------------------

Each benchmark typically reports:

* **time** – execution time in seconds
* **ratio** – comparison factor between two commits
* **std** – standard deviation across runs
* **samples** – count of timing samples

General guidelines:

* Treat **ratio ≥ 1.2** as a likely regression.
* Large **std** indicates noise; rerun with ``asv run``.
* For unstable benchmarks, consider running more samples.

------------------------------
Viewing Results with ASV Preview
------------------------------

To launch an interactive viewer:

.. code-block:: bash

    asv preview

This opens a local web interface showing timing graphs, regressions,
benchmark descriptions, and environment details.

-----------------------------------------------
Detecting and Reporting Regressions in Pull Requests
-----------------------------------------------

Before submitting a PR, contributors should check for regressions.

Recommended process:

1. Build SymPy in-place:

   .. code-block:: bash

       python setup.py build_ext --inplace

2. Run relevant benchmarks:

   .. code-block:: bash

       asv run benchmarks/<module>

3. Compare with ``master``:

   .. code-block:: bash

       asv compare master HEAD

If a regression is found:

* Determine whether it is expected.
* If unexpected, profile or isolate the change.
* Mention ASV results in the PR description.

----------------------
Useful References
----------------------

* ASV documentation: https://asv.readthedocs.io/en/stable/
* SymPy benchmark suite: ``benchmarks/`` directory
* Contributor Guide: :ref:`contributing`
