#!/bin/bash
python3 -m pytest --strict -ra --veryquickcheck --pep8 --flakes --durations 10
if [[ $? -ne 0 ]]; then
    python3 -m pytest --strict --lf --verbose
else
    python3 -m pytest --strict -ra --quickcheck --durations 10
fi
