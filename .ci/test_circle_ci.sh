#!/bin/bash
TESTCMD=python3 -m pytest --strict -ra
for fname in $(git diff --name-only master); do
    if [[ $fname == *.py ]]; then
        $TESTCMD --doctest-modules $fname
    elif [[ $fname == *.rst ]]; then
        python3 -m doctest $fname
    fi
done
$TESTCMD --veryquickcheck --pep8 --flakes --durations 10 >veryquick.log
tail -n 80 veryquick.log
if [[ $? -ne 0 ]]; then
    python3 -m pytest --strict --lf --verbose
else
    python3 -m pytest --strict -ra --quickcheck --durations 10
fi
