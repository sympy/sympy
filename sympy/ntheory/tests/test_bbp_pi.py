from sympy.ntheory.bbp_pi import pi_hex_digits


def test_hex_pi_nth_digits():
    assert pi_hex_digits(0) == '3243f6a8885a30'
    assert pi_hex_digits(1) == '243f6a8885a308'
    assert pi_hex_digits(10000) == '68ac8fcfb8016c'
