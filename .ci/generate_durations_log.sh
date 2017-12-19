#!/bin/bash
ABS_REPO_PATH=$(unset CDPATH && cd "$(dirname "$0")/.." && echo $PWD)
python3 -m pytest --durations 0 >$ABS_REPO_PATH/.ci/durations.log
