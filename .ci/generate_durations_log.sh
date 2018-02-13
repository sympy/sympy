#!/bin/bash
ABS_REPO_PATH=$(unset CDPATH && cd "$(dirname "$0")/.." && echo $PWD)
python3 -m pytest -ra -k "not test_bicycle" --durations 0 >$ABS_REPO_PATH/.ci/durations.log
cat <<EOF >${ABS_REPO_PATH}/blacklisted.json
{
    "sympy/physics/mechanics/tests/test_kane3.py": [
        "test_bicycle"
    ],
}
EOF
