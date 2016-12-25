import sys

from sympy.interactive.printing import _init_python_printing


def test_init_python_printing():
    original_display_hook = sys.displayhook

    def new_display_hook(*_):
        new_display_hook.called = True

    sys.displayhook = new_display_hook

    try:
        _init_python_printing(str)
        # check for correct overwriting of displayhook
        assert getattr(sys.displayhook, 'is_sympy_displayhook', False)
        # call displayhook with none-sympy object
        sys.displayhook(object())
        # check whether original display hook has been called
        assert getattr(new_display_hook, 'called', False)
    finally:
        sys.displayhook = original_display_hook
