===============================
Performance Benchmarking Guide
===============================

This guide explains how to use **Airspeed Velocity (asv)** for benchmarking SymPy’s performance.
It is useful for contributors who want to measure runtime changes, compare commits, or detect regressions.

-------------------------------------
1. Clone the SymPy benchmarks repository
-------------------------------------

SymPy’s performance benchmarks are maintained separately:

    https://github.com/sympy/sympy_benchmarks

Clone it alongside your local SymPy repository:

.. code-block:: bash

    git clone https://github.com/sympy/sympy_benchmarks.git
    cd sympy_benchmarks

-------------------------------------
2. Install dependencies
-------------------------------------

Install **asv**:

.. code-block:: bash

    pip install asv

Optionally, install your development version of SymPy:

.. code-block:: bash

    pip install -e ../sympy

-------------------------------------
3. Run benchmarks
-------------------------------------

Run all available benchmarks:

.. code-block:: bash

    asv run

Run only a specific group (for example, integration benchmarks):

.. code-block:: bash

    asv run --bench "integrate"

Compare two branches or commits:

.. code-block:: bash

    asv compare master HEAD

-------------------------------------
4. View results
-------------------------------------

To generate and open a local benchmark dashboard:

.. code-block:: bash

    asv publish
    asv preview

-------------------------------------
5. Interpret results
-------------------------------------

* **Mean** – Average time per benchmark run  
* **Std** – Standard deviation (timing variation)  
* **Ratio** – Performance change relative to a baseline

A ratio > 1 means slower performance; < 1 means faster.

-------------------------------------
6. Reporting regressions
-------------------------------------

If you observe a slowdown:
1. Verify using `asv run` on both commits.  
2. Record the commit hashes and benchmark names.  
3. Report in an issue with summary data.

-------------------------------------
7. Further reading
-------------------------------------

* SymPy Benchmarks Repository: https://github.com/sympy/sympy_benchmarks
* Airspeed Velocity Docs: https://asv.readthedocs.io/
