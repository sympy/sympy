def test_rounding():
    assert S(123).round(-2) == 100
    assert S(123.456).round(2) == 123.46
