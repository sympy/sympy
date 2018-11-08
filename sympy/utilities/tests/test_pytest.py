from sympy.utilities.pytest import raises, USE_PYTEST

if USE_PYTEST:
    import py.test
    pytestmark = py.test.mark.skipif(USE_PYTEST,
                                     reason=("using py.test"))

# Test callables


def test_expected_exception_is_silent_callable():
    def f():
        raise ValueError()
    raises(ValueError, f)


def test_lack_of_exception_triggers_AssertionError_callable():
    try:
        raises(Exception, lambda: 1 + 1)
        assert False
    except AssertionError as e:
        assert str(e) == "DID NOT RAISE"


def test_unexpected_exception_is_passed_through_callable():
    def f():
        raise ValueError("some error message")
    try:
        raises(TypeError, f)
        assert False
    except ValueError as e:
        assert str(e) == "some error message"

# Test with statement

def test_expected_exception_is_silent_with():
    with raises(ValueError):
        raise ValueError()


def test_lack_of_exception_triggers_AssertionError_with():
    try:
        with raises(Exception):
            1 + 1
        assert False
    except AssertionError as e:
        assert str(e) == "DID NOT RAISE"


def test_unexpected_exception_is_passed_through_with():
    try:
        with raises(TypeError):
            raise ValueError("some error message")
        assert False
    except ValueError as e:
        assert str(e) == "some error message"

# Now we can use raises() instead of try/catch
# to test that a specific exception class is raised


def test_second_argument_should_be_callable_or_string():
    raises(TypeError, lambda: raises("irrelevant", 42))
