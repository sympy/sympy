#!/bin/bash -e
ABS_REPO_PATH=$(unset CDPATH && cd "$(dirname "$0")/.." && echo $PWD)
cat <<EOF >${ABS_REPO_PATH}/.ci/blacklisted.json
{
    "sympy/utilities/tests/test_wester.py": [
        "test_W25"
    ]
}
EOF
${PYTHON:-python} -m pytest -ra --durations 0 --verbose | tee $ABS_REPO_PATH/.ci/durations.log
