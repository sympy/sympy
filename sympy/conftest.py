import sys
sys._running_pytest = True

import pytest
from sympy.core.cache import clear_cache

def pytest_report_header(config):
    from sympy.utilities.misc import ARCH
    s = "architecture: %s\n" % ARCH
    from sympy.core.cache import USE_CACHE
    s += "cache:        %s\n" % USE_CACHE
    from sympy.polys.domains import GROUND_TYPES
    s += "ground types: %s\n" % GROUND_TYPES
    return s

def pytest_addoption(parser):
    parser.addoption("--slow", dest="runslow", action="store_true",
        help="allow slow tests to run")

def pytest_configure(config):
    # register an additional marker
    config.addinivalue_line("markers", "slow: slow test")

def pytest_runtest_setup(item):
    if not isinstance(item, pytest.Function):
        return
    if not item.config.getvalue("runslow") and hasattr(item.obj, 'slow'):
        pytest.skip("slow test: pass --slow to run")


def pytest_terminal_summary(terminalreporter):
    if (terminalreporter.stats.get('error', None) or
            terminalreporter.stats.get('failed', None)):
        terminalreporter.write_sep(' ', 'DO *NOT* COMMIT!', red=True, bold=True)

def pytest_runtest_teardown():
    clear_cache()
