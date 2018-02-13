from sympy.utilities.pytest import raises, USE_PYTEST, assert_raise_message

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


def test_assert_raise_message():
    import re

    message1 = re.compile(r"division by zero", re.IGNORECASE)
    message2 = re.compile(r"zero", re.IGNORECASE)
    message3 = re.compile(r"one", re.IGNORECASE)
    message4 = re.compile(r"some error message", re.IGNORECASE)

    assert_raise_message(ZeroDivisionError, message1, lambda: 1/0) # error message: 'division by zero'
    assert_raise_message(ZeroDivisionError, message2, lambda: 1/0)

    raises(AssertionError, lambda: assert_raise_message(ZeroDivisionError, message3, lambda: 1/0))

    #with assert_raise_message(ValueError, message4) as exc_wrapper:
    #    raise ValueError("some error message")
    #assert isinstance(exc_wrapper, ValueError)
