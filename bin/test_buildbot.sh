#!/usr/bin/env bash

# We are using the NiPy buildbot to run our time intensive tests: with all
# optional dependencies, without caching, and tests marked as slow.

# Exit on error
set -e
# Echo each command
set -x

# run all tests (all install dependencies should be installed)
python - <<'EOF'
import sympy
if not sympy.test():
    raise Exception('Tests failed')
EOF
# run only the slow tests, time out after 30 minutes
python - <<'EOF'
import sympy
if not sympy.test(slow=True, timeout=1800):
    raise Exception('Tests failed')
EOF

# run all the tests without the cache
export SYMPY_USE_CACHE=no
python - <<'EOF'
import sympy
if not sympy.test():
    raise Exception('Tests failed')
EOF
python - <<'EOF'
import sympy
if not sympy.test(slow=True, timeout=1800):
    raise Exception('Tests failed')
EOF
