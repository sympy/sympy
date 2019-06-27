This folder currently holds ``durations.json`` and ``blacklisted.json`` which are used in our pytest configuration.
``durations.json`` allows the two options ``--quickcheck`` and ``--veryquickcheck`` to be passed to ``pytest``.
Currently, the list of slow tests is updated intermittently by running::

  $ ./.ci/generate_durations_log.sh
  $ ./.ci/parse_durations_log.py
  $ git commit -am "pytest: Updated .ci/durations.json"

make sure you have a C compiler, a Fortran compiler, numpy, scipy, cython, matplotlib & sage installed (the latter can
be tricky, in which case just make sure that no sage tests are removed in the git diff).
