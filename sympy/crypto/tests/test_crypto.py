from sympy.crypto.crypto import (alphabet_of_cipher, cycle_list,
      encipher_shift, encipher_affine, encipher_substitution,
      encipher_vigenere, decipher_vigenere, matrix_inverse_mod,
      encipher_hill, decipher_hill, encipher_bifid5, encipher_bifid6,
      bifid5_square, bifid6_square, bifid7_square,
      encipher_bifid7, decipher_bifid5, decipher_bifid6, encipher_kid_rsa,
      decipher_kid_rsa, kid_rsa_private_key, kid_rsa_public_key,
      decipher_rsa, rsa_private_key, rsa_public_key, encipher_rsa,
      lfsr_connection_polynomial, lfsr_autocorrelation, lfsr_sequence)
from sympy import Matrix
from sympy.polys.domains.pythonfinitefield import PythonFiniteField as FF
from sympy import symbols

def test_alphabet_of_cipher():
    assert alphabet_of_cipher()[0] == "A"
    assert alphabet_of_cipher(symbols = "1z") == ["1", "z"]


def test_cycle_list():
    assert cycle_list(3,4) == [3, 0, 1, 2]
    assert cycle_list(-1,4) == [3, 0, 1, 2]
    assert cycle_list(1,4) == [1, 2, 3, 0]

def test_encipher_shift():
    assert encipher_shift("ABC", 0)=="ABC"
    assert encipher_shift("ABC", 1)=="BCD"
    assert encipher_shift("ABC", -1)=="ZAB"


def test_encipher_affine():
    assert encipher_affine("ABC", (1, 0))=="ABC"
    assert encipher_affine("ABC", (1, 1))=="BCD"
    assert encipher_affine("ABC", (-1, 0))=="AZY"
    assert encipher_affine("ABC", (-1, 1), symbols = "ABCD")=="BAD"
    assert encipher_affine("123", (-1, 1), symbols = "1234")=="214"


def test_encipher_substitution():
    assert encipher_substitution("ABC", "BAC", symbols = "ABC")=="BAC"
    assert encipher_substitution("123", "124", symbols = "1234")=="124"


def test_encipher_vigenere():
    assert encipher_vigenere("ABC", "ABC")=="ACE"
    assert encipher_vigenere("ABC", "ABC", symbols = "ABCD")=="ACA"
    assert encipher_vigenere("ABC", "AB", symbols = "ABCD")=="ACC"
    assert encipher_vigenere("AB", "ABC", symbols = "ABCD")=="AC"
    assert encipher_vigenere("A", "ABC", symbols = "ABCD")=="A"

def test_decipher_vigenere():
    assert decipher_vigenere("ABC", "ABC")=="AAA"
    assert decipher_vigenere("ABC", "ABC", symbols = "ABCD")=="AAA"
    assert decipher_vigenere("ABC", "AB", symbols = "ABCD")=="AAC"
    assert decipher_vigenere("AB", "ABC", symbols = "ABCD")=="AA"
    assert decipher_vigenere("A", "ABC", symbols = "ABCD")=="A"


def test_matrix_inverse_mod():
    A = Matrix(2, 2, [1, 2, 3, 4])
    Ai = Matrix(2, 2, [1, 1, 0, 1])
    assert matrix_inverse_mod(A, 3)==Ai
    A = Matrix(2, 2, [1, 0, 0, 1])
    assert matrix_inverse_mod(A, 2)==A


def test_encipher_hill():
    A = Matrix(2, 2, [1, 2, 3, 5])
    assert encipher_hill("ABCD", A)=="CFIV"
    A = Matrix(2, 2, [1, 0, 0, 1])
    assert encipher_hill("ABCD", A)=="ABCD"
    assert encipher_hill("ABCD", A, symbols = "ABCD")=="ABCD"
    A = Matrix(2, 2, [1, 2, 3, 5])
    assert encipher_hill("ABCD", A, symbols="ABCD")=="CBAB"
    assert encipher_hill("AB", A, symbols = "ABCD")=="CB"


def test_decipher_hill():
    A = Matrix(2, 2, [1, 2, 3, 5])
    assert decipher_hill("CFIV", A)=="ABCD"
    A = Matrix(2, 2, [1, 0, 0, 1])
    assert decipher_hill("ABCD", A)=="ABCD"
    assert decipher_hill("ABCD", A, symbols = "ABCD")=="ABCD"
    A = Matrix(2, 2, [1, 2, 3, 5])
    assert decipher_hill("CBAB", A, symbols="ABCD")=="ABCD"
    assert decipher_hill("CB", A, symbols = "ABCD")=="AB"


def test_encipher_bifid5():
    assert encipher_bifid5("AB", "AB")=="AB"
    assert encipher_bifid5("AB", "CD")=="CO"
    assert encipher_bifid5("ab", "c")=="CH"
    assert encipher_bifid5("a bc", "b")=="BAC"


def test_bifid5_square():
    A = alphabet_of_cipher()
    A.remove("J")
    f = lambda i,j: symbols(A[5*i+j])
    M = Matrix(5, 5, f)
    assert bifid5_square("")==M


def test_decipher_bifid5():
    assert decipher_bifid5("AB", "AB")=="AB"
    assert decipher_bifid5("CO", "CD")=="AB"
    assert decipher_bifid5("ch", "c")=="AB"
    assert decipher_bifid5("b ac", "b")=="ABC"


def test_bifid7_square():
    A = alphabet_of_cipher()+[str(a) for a in range(23)]
    f = lambda i,j: symbols(A[7*i+j])
    M = Matrix(7, 7, f)
    assert bifid7_square("")==M


def test_encipher_bifid7():
    assert encipher_bifid7("AB", "AB")=="AB"
    assert encipher_bifid7("AB", "CD")=="CR"
    assert encipher_bifid7("ab", "c")=="CJ"
    assert encipher_bifid7("a bc", "b")=="BAC"


def test_encipher_bifid6():
    assert encipher_bifid6("AB", "AB")=="AB"
    assert encipher_bifid6("AB", "CD")=="CP"
    assert encipher_bifid6("ab", "c")=="CI"
    assert encipher_bifid6("a bc", "b")=="BAC"


def test_decipher_bifid6():
    assert decipher_bifid6("AB", "AB")=="AB"
    assert decipher_bifid6("CP", "CD")=="AB"
    assert decipher_bifid6("ci", "c")=="AB"
    assert decipher_bifid6("b ac", "b")=="ABC"


def test_bifid6_square():
    A = alphabet_of_cipher()+[str(a) for a in range(10)]
    f = lambda i,j: symbols(A[6*i+j])
    M = Matrix(6, 6, f)
    assert bifid6_square("")==M


def test_rsa_public_key():
    assert rsa_public_key(2,2,1)==(4,1)
    assert rsa_public_key(2,3,1)==(6,1)
    assert rsa_public_key(5,3,3)==(15,3)


def test_rsa_private_key():
    assert rsa_private_key(2,2,1)==(4,1)
    assert rsa_private_key(2,3,1)==(6,1)
    assert rsa_private_key(5,3,3)==(15,3)


def test_encipher_rsa():
    assert encipher_rsa(2,2,1,2)==2
    assert encipher_rsa(2,3,1,2)==2
    assert encipher_rsa(5,3,3,2)==8


def test_decipher_rsa():
    assert decipher_rsa(2,2,1,2)==2
    assert decipher_rsa(2,3,1,2)==2
    assert decipher_rsa(5,3,3,8)==2


def test_kid_rsa_public_key():
    assert kid_rsa_public_key(1,2,1,1)==(5,2)
    assert kid_rsa_public_key(1,2,2,1)==(8,3)
    assert kid_rsa_public_key(1,2,1,2)==(7,2)


def test_kid_rsa_private_key():
    assert kid_rsa_private_key(1,2,1,1)==(5,3)
    assert kid_rsa_private_key(1,2,2,1)==(8,3)
    assert kid_rsa_private_key(1,2,1,2)==(7,4)


def test_encipher_kid_rsa():
    assert encipher_kid_rsa(1,(5,2))==2
    assert encipher_kid_rsa(1,(8,3))==3
    assert encipher_kid_rsa(1,(7,2))==2


def test_decipher_kid_rsa():
    assert decipher_kid_rsa(2,(5,3))==1
    assert decipher_kid_rsa(3,(8,3))==1
    assert decipher_kid_rsa(2,(7,4))==1


def test_lfsr_sequence():
    F = FF(2)
    assert lfsr_sequence([F(1)], [F(1)], 2)==[F(1), F(1)]
    assert lfsr_sequence([F(0)], [F(1)], 2)==[F(1), F(0)]
    F = FF(3)
    assert lfsr_sequence([F(1)], [F(1)], 2)==[F(1), F(1)]
    assert lfsr_sequence([F(0)], [F(2)], 2)==[F(2), F(0)]
    assert lfsr_sequence([F(1)], [F(2)], 2)==[F(2), F(2)]


def test_lfsr_autocorrelation():
    F = FF(2)
    s = lfsr_sequence([F(1),F(0)], [F(0),F(1)], 5)
    assert lfsr_autocorrelation(s,2,0)==1
    assert lfsr_autocorrelation(s,2,1)==-1


def test_lfsr_connection_polynomial():
    F = FF(2)
    x = symbols("x")
    s = lfsr_sequence([F(1),F(0)], [F(0),F(1)], 5)
    assert lfsr_connection_polynomial(s)==x**2+1
