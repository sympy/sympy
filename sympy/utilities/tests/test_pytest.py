from sympy.utilities.pytest import raises

# Test callables

def test_expected_exception_is_silent():
    def f():
        raise ValueError()
    raises(ValueError, f)

def test_lack_of_exception_triggers_AssertionError():
    try:
        raises(Exception, lambda: 1+1)
        assert False
    except AssertionError, e:
        assert str(e) == "DID NOT RAISE"

def test_unexpected_exception_is_passed_through():
    def f():
        raise ValueError("some error message")
    try:
        raises(TypeError, f)
        assert False
    except ValueError, e:
        assert str(e) == "some error message"

# Test compilable strings

def test_expected_exception_is_silent():
    raises(ValueError, "raise ValueError()")

def test_lack_of_exception_triggers_AssertionError():
    try:
        raises(Exception, "1+1")
        assert False
    except AssertionError, e:
        assert str(e) == "DID NOT RAISE"

def test_unexpected_exception_is_passed_through():
    try:
        raises(TypeError, 'raise ValueError("some error message")')
        assert False
    except ValueError, e:
        assert str(e) == "some error message"

# Test with statement

def test_expected_exception_is_silent():
    with raises(ValueError):
        raise ValueError()

def test_lack_of_exception_triggers_AssertionError():
    try:
        with raises(Exception):
            1+1
        assert False
    except AssertionError, e:
        assert str(e) == "DID NOT RAISE"

def test_unexpected_exception_is_passed_through():
    try:
        with raises(TypeError):
            raise ValueError("some error message")
        assert False
    except ValueError, e:
        assert str(e) == "some error message"

# Now we can use raises() instead of try/catch
# to test that a specific exception class is raised

def test_second_argument_should_be_callable_or_string():
    raises(TypeError, lambda: raises("irrelevant", 42))
