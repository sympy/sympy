#!/bin/bash
REPOROOT=$(unset CDPATH && cd "$(dirname "$0")/.." && echo $PWD)
for f in $(find $REPOROOT/sympy -name "*.py"); do
    output=$(python3 -We:invalid -m compileall -f $f)
    if [ $? -ne 0 ]; then
       echo $output
       break;
    fi
done
