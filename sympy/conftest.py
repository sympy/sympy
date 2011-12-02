import sys
sys._running_pytest = True

from sympy.core.cache import clear_cache

def pytest_report_header(config):
    from sympy.utilities.misc import ARCH
    s = "architecture: %s\n" % ARCH
    from sympy.core.cache import USE_CACHE
    s += "cache:        %s\n" % USE_CACHE
    from sympy.polys.domains import GROUND_TYPES
    s += "ground types: %s\n" % GROUND_TYPES
    return s

def pytest_terminal_summary(terminalreporter):
    if (terminalreporter.stats.get('error', None) or
            terminalreporter.stats.get('failed', None)):
        terminalreporter.write_sep(' ', 'DO *NOT* COMMIT!', red=True, bold=True)

def pytest_runtest_teardown():
    clear_cache()
