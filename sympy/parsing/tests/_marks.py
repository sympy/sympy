import subprocess
import pytest


def _skip_if_callable(args):
    """ A mark to use to skip a test if a command is callable
    """
    try:
        subprocess.Popen(args, shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE).communicate()
        has_cmd = True
    except FileNotFoundError:
        has_cmd = False

    return [
        pytest.mark.skipif(has_cmd,
                           reason="COULD call '%s'".format(" ".join(args))),
        pytest.mark.skipif(not has_cmd,
                           reason="COULDN'T call '%s'".format(" ".join(args)))
    ]


def _skip_if_importable(module_name):
    try:
        __import__(module_name)
        has_module = True
    except ImportError:
        has_module = False

    return [
        pytest.mark.skipif(has_module,
                           reason="COULD import '%s'".format(module_name)),
        pytest.mark.skipif(not has_module,
                           reason="COULDN'T import '%s'".format(module_name))
    ]

skip_if_has_antlr4, skip_if_missing_antlr4 = _skip_if_callable(["antlr4"])

skip_if_has_antlr4runtime, skip_if_missing_antlr4runtime = _skip_if_importable("antlr4")
