#!/bin/bash -x
TESTCMD="python3 -m pytest --strict"
for fname in $(git diff --name-only master); do
    if [[ $fname == *.py ]]; then
        $TESTCMD --doctest-modules $fname
    elif [[ $fname == *.rst ]]; then
        python3 -m doctest $fname
    fi
done
$TESTCMD -ra --veryquickcheck --pep8 --flakes --durations 10 >veryquick.log
if [[ $? -ne 0 ]]; then
    $TESTCMD --lf --verbose
else
    tail -n 80 veryquick.log
fi
