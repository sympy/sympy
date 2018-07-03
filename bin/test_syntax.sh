#!/bin/bash
REPOROOT=$(unset CDPATH && cd "$(dirname "$0")/.." && echo $PWD)
for f in $(find $REPOROOT/sympy -name "*.py"); do
    output=$(python3 -We:invalid -m py_compile $f)
    if [ $? -ne 0 ]; then
       >&2 echo $output
       break;
    fi
done
