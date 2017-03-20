# -*- coding: utf-8 -*-

from __future__ import print_function, division, absolute_import

import os
import json

import pytest



durations_path = os.path.join('.ci', 'durations.json')
if os.path.exists(durations_path):
    veryslow_group, slow_group = json.loads(open(durations_path, 'rt').read())
else:
    veryslow_group, slow_group = None, None


def pytest_addoption(parser):
    parser.addoption("--slow", dest="runslow", action="store_true",
                     help="allow slow tests to run")
    parser.addoption("--veryslow", dest="runveryslow", action="store_true",
                     help="allow very slow tests to run")


def pytest_configure(config):
    # register an additional marker
    config.addinivalue_line("markers", "slow: slow test")
    config.addinivalue_line("markers", "veryslow: very slow test")


def pytest_runtest_setup(item):
    if not isinstance(item, pytest.Function):
        return

    if veryslow_group is not None:
        if item.nodeid in veryslow_group:
            pytest.skip("very slow test: pass --veryslow to run")
            return

    if slow_group is not None:
        if item.nodeid in slow_group:
            pytest.skip("slow test: pass --slow to run")
            return

    if not item.config.getvalue("runslow") and hasattr(item.obj, 'slow'):
        pytest.skip("slow test: pass --slow to run")
        return
    if not item.config.getvalue("runveryslow") and hasattr(item.obj, 'veryslow'):
        pytest.skip("very slow test: pass --veryslow to run")
        return
