import sys
sys._running_pytest = True

from sympy.core.cache import clear_cache

def pytest_terminal_summary(terminalreporter):
    if (terminalreporter.stats.get('error', None) or
            terminalreporter.stats.get('failed', None)):
        terminalreporter.write_sep(' ', 'DO *NOT* COMMIT!', red=True, bold=True)

def pytest_runtest_teardown():
    clear_cache()
