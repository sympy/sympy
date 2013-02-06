#! /usr/bin/env bash

set -e

if [[ "${TEST_RST}" == "true" ]]; then

bin/doctest

else

# We change directories to make sure that python won't find the copy
# of sympy in the source directory.
mkdir empty
cd empty
cat << EOF | python
import sympy
t1=sympy.test()
t2=sympy.doctest()
if not (t1 and t2):
    raise Exception('Tests failed')
EOF

fi
