.. _benchmarking-guide:

Performance Benchmarking with ASV
=================================

SymPy uses `Airspeed Velocity (ASV) <https://asv.readthedocs.io/>`_ to measure
performance and detect regressions over time. This guide explains how
contributors can run benchmarks locally, compare results across commits,
interpret ASV output, and use SymPy’s existing benchmark repositories.

Why Benchmark?
--------------

Small code changes may significantly affect performance in symbolic
manipulation, simplification, assumptions, or expression rewriting.
Benchmarking helps contributors:

* measure the performance impact of changes,
* identify regressions before submitting a pull request,
* justify optimizations with quantitative data.

ASV Installation
----------------

Install ASV into any environment:

.. code-block:: bash

    pip install asv

Verify installation:

.. code-block:: bash

    asv --help


Benchmark Locations
-------------------

SymPy benchmarks exist in two places:

1. The ``benchmarks/`` directory inside the main SymPy repository.
2. The external repository: https://github.com/sympy/sympy_benchmarks

The external repository contains a larger, curated benchmark suite and
historical performance results. It is useful for deeper investigations and
for comparing local results with published ASV histories.


Running Benchmarks
------------------

Run all benchmarks:

.. code-block:: bash

    asv run

Run a specific benchmark file:

.. code-block:: bash

    asv run benchmarks/bench_assumptions.py

Run a specific benchmark class or method:

.. code-block:: bash

    asv run "assumptions.AssumptionBench.time_simplify"

Run a quick subset (useful during development):

.. code-block:: bash

    asv run --quick


Comparing Performance Across Commits
------------------------------------

This is the most important ASV command for contributors.

Compare performance between ``master`` and the current branch:

.. code-block:: bash

    asv compare master HEAD

Compare arbitrary commits:

.. code-block:: bash

    asv compare <commit1> <commit2>

Interpreting ASV Output
-----------------------

Example comparison:

::

    benchmark                                ratio    ci   std
    -------------------------------------    ------   ---  ----
    bench_simplify.time_simplify            1.23x    5%   0.01

Meaning:

* **ratio > 1** → slower (possible regression)
* **ratio < 1** → faster (improvement)
* **ci** = confidence interval
* **std** = run variability

A slowdown above ~5% should be investigated and mentioned in the pull request.


Previewing Benchmark Reports
----------------------------

Generate a local HTML performance report:

.. code-block:: bash

    asv preview

This opens an interactive ASV dashboard showing performance history and
benchmark details.


Detecting Regressions Before Submitting a PR
--------------------------------------------

Recommended workflow:

1. Benchmark ``master``:

   .. code-block:: bash

       git checkout master
       asv run

2. Benchmark your branch:

   .. code-block:: bash

       git checkout my-branch
       asv run

3. Compare:

   .. code-block:: bash

       asv compare master my-branch

4. Include a short summary in your PR:

   Example::

       ASV comparison: no regressions observed.
       Largest deviation: bench_simplify.time_simplify = 0.97x (3% faster).


Tips for Reliable Benchmarking
------------------------------

* Close CPU-intensive programs.
* Run benchmarks multiple times if results appear noisy.
* Use ``asv run --quick`` during development.
* Use full ASV runs before finalizing a PR.


Additional Resources
--------------------

**sympy_benchmarks repository**  
A larger curated benchmark suite with historical results:

* https://github.com/sympy/sympy_benchmarks

Contains:

* detailed README with setup instructions,
* examples of high-quality benchmarks,
* long-term performance history.

Contributors are encouraged to consult it when analyzing tricky performance
issues.

**ASV Documentation**  
Comprehensive ASV user guide:

* https://asv.readthedocs.io/en/stable/
Covers advanced topics like custom configurations, CI integration, and
report generation.