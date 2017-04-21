from sympy.core import symbols
from sympy.core.compatibility import range
from sympy.crypto.crypto import (cycle_list,
      encipher_shift, encipher_affine, encipher_substitution,
      check_and_join, encipher_vigenere, decipher_vigenere,
      encipher_hill, decipher_hill, encipher_bifid5, encipher_bifid6,
      bifid5_square, bifid6_square, bifid5, bifid6, bifid10,
      decipher_bifid5, decipher_bifid6, encipher_kid_rsa,
      decipher_kid_rsa, kid_rsa_private_key, kid_rsa_public_key,
      decipher_rsa, rsa_private_key, rsa_public_key, encipher_rsa,
      lfsr_connection_polynomial, lfsr_autocorrelation, lfsr_sequence,
      encode_morse, decode_morse, elgamal_private_key, elgamal_public_key,
      encipher_elgamal, decipher_elgamal, dh_private_key, dh_public_key,
      dh_shared_key, decipher_shift, decipher_affine, encipher_bifid,
      decipher_bifid, bifid_square, padded_key, uniq, playfair_matrix,
      encipher_playfair, decipher_playfair, mix_columns, expand_key,
      asciify, hexify, binify, encipher_rijndael,decipher_rijndael)
from sympy.matrices import Matrix
from sympy.ntheory import isprime, is_primitive_root
from sympy.polys.domains import FF

from sympy.utilities.pytest import raises, slow

from random import randrange


def test_cycle_list():
    assert cycle_list(3, 4) == [3, 0, 1, 2]
    assert cycle_list(-1, 4) == [3, 0, 1, 2]
    assert cycle_list(1, 4) == [1, 2, 3, 0]


def test_encipher_shift():
    assert encipher_shift("ABC", 0) == "ABC"
    assert encipher_shift("ABC", 1) == "BCD"
    assert encipher_shift("ABC", -1) == "ZAB"
    assert decipher_shift("ZAB", -1) == "ABC"


def test_encipher_affine():
    assert encipher_affine("ABC", (1, 0)) == "ABC"
    assert encipher_affine("ABC", (1, 1)) == "BCD"
    assert encipher_affine("ABC", (-1, 0)) == "AZY"
    assert encipher_affine("ABC", (-1, 1), symbols="ABCD") == "BAD"
    assert encipher_affine("123", (-1, 1), symbols="1234") == "214"
    assert encipher_affine("ABC", (3, 16)) == "QTW"
    assert decipher_affine("QTW", (3, 16)) == "ABC"


def test_encipher_substitution():
    assert encipher_substitution("ABC", "BAC", "ABC") == "BAC"
    assert encipher_substitution("123", "1243", "1234") == "124"


def test_check_and_join():
    assert check_and_join("abc") == "abc"
    assert check_and_join(uniq("aaabc")) == "abc"
    assert check_and_join("ab c".split()) == "abc"
    assert check_and_join("abc", "a", filter=True) == "a"
    raises(ValueError, lambda: check_and_join('ab', 'a'))


def test_encipher_vigenere():
    assert encipher_vigenere("ABC", "ABC") == "ACE"
    assert encipher_vigenere("ABC", "ABC", symbols="ABCD") == "ACA"
    assert encipher_vigenere("ABC", "AB", symbols="ABCD") == "ACC"
    assert encipher_vigenere("AB", "ABC", symbols="ABCD") == "AC"
    assert encipher_vigenere("A", "ABC", symbols="ABCD") == "A"


def test_decipher_vigenere():
    assert decipher_vigenere("ABC", "ABC") == "AAA"
    assert decipher_vigenere("ABC", "ABC", symbols="ABCD") == "AAA"
    assert decipher_vigenere("ABC", "AB", symbols="ABCD") == "AAC"
    assert decipher_vigenere("AB", "ABC", symbols="ABCD") == "AA"
    assert decipher_vigenere("A", "ABC", symbols="ABCD") == "A"


def test_encipher_hill():
    A = Matrix(2, 2, [1, 2, 3, 5])
    assert encipher_hill("ABCD", A) == "CFIV"
    A = Matrix(2, 2, [1, 0, 0, 1])
    assert encipher_hill("ABCD", A) == "ABCD"
    assert encipher_hill("ABCD", A, symbols="ABCD") == "ABCD"
    A = Matrix(2, 2, [1, 2, 3, 5])
    assert encipher_hill("ABCD", A, symbols="ABCD") == "CBAB"
    assert encipher_hill("AB", A, symbols="ABCD") == "CB"
    # message length, n, does not need to be a multiple of k;
    # it is padded
    assert encipher_hill("ABA", A) == "CFGC"
    assert encipher_hill("ABA", A, pad="Z") == "CFYV"


def test_decipher_hill():
    A = Matrix(2, 2, [1, 2, 3, 5])
    assert decipher_hill("CFIV", A) == "ABCD"
    A = Matrix(2, 2, [1, 0, 0, 1])
    assert decipher_hill("ABCD", A) == "ABCD"
    assert decipher_hill("ABCD", A, symbols="ABCD") == "ABCD"
    A = Matrix(2, 2, [1, 2, 3, 5])
    assert decipher_hill("CBAB", A, symbols="ABCD") == "ABCD"
    assert decipher_hill("CB", A, symbols="ABCD") == "AB"
    # n does not need to be a multiple of k
    assert decipher_hill("CFA", A) == "ABAA"


def test_encipher_bifid5():
    assert encipher_bifid5("AB", "AB") == "AB"
    assert encipher_bifid5("AB", "CD") == "CO"
    assert encipher_bifid5("ab", "c") == "CH"
    assert encipher_bifid5("a bc", "b") == "BAC"


def test_bifid5_square():
    A = bifid5
    f = lambda i, j: symbols(A[5*i + j])
    M = Matrix(5, 5, f)
    assert bifid5_square("") == M


def test_decipher_bifid5():
    assert decipher_bifid5("AB", "AB") == "AB"
    assert decipher_bifid5("CO", "CD") == "AB"
    assert decipher_bifid5("ch", "c") == "AB"
    assert decipher_bifid5("b ac", "b") == "ABC"


def test_encipher_bifid6():
    assert encipher_bifid6("AB", "AB") == "AB"
    assert encipher_bifid6("AB", "CD") == "CP"
    assert encipher_bifid6("ab", "c") == "CI"
    assert encipher_bifid6("a bc", "b") == "BAC"


def test_decipher_bifid6():
    assert decipher_bifid6("AB", "AB") == "AB"
    assert decipher_bifid6("CP", "CD") == "AB"
    assert decipher_bifid6("ci", "c") == "AB"
    assert decipher_bifid6("b ac", "b") == "ABC"


def test_bifid6_square():
    A = bifid6
    f = lambda i, j: symbols(A[6*i + j])
    M = Matrix(6, 6, f)
    assert bifid6_square("") == M


def test_rsa_public_key():
    assert rsa_public_key(2, 2, 1) == (4, 1)
    assert rsa_public_key(2, 3, 1) == (6, 1)
    assert rsa_public_key(5, 3, 3) == (15, 3)
    assert rsa_public_key(8, 8, 8) is False


def test_rsa_private_key():
    assert rsa_private_key(2, 2, 1) == (4, 1)
    assert rsa_private_key(2, 3, 1) == (6, 1)
    assert rsa_private_key(5, 3, 3) == (15, 3)
    assert rsa_private_key(23,29,5) == (667,493)
    assert rsa_private_key(8, 8, 8) is False


def test_encipher_rsa():
    puk = rsa_public_key(2, 2, 1)
    assert encipher_rsa(2, puk) == 2
    puk = rsa_public_key(2, 3, 1)
    assert encipher_rsa(2, puk) == 2
    puk = rsa_public_key(5, 3, 3)
    assert encipher_rsa(2, puk) == 8


def test_decipher_rsa():
    prk = rsa_private_key(2, 2, 1)
    assert decipher_rsa(2, prk) == 2
    prk = rsa_private_key(2, 3, 1)
    assert decipher_rsa(2, prk) == 2
    prk = rsa_private_key(5, 3, 3)
    assert decipher_rsa(8, prk) == 2


def test_kid_rsa_public_key():
    assert kid_rsa_public_key(1, 2, 1, 1) == (5, 2)
    assert kid_rsa_public_key(1, 2, 2, 1) == (8, 3)
    assert kid_rsa_public_key(1, 2, 1, 2) == (7, 2)


def test_kid_rsa_private_key():
    assert kid_rsa_private_key(1, 2, 1, 1) == (5, 3)
    assert kid_rsa_private_key(1, 2, 2, 1) == (8, 3)
    assert kid_rsa_private_key(1, 2, 1, 2) == (7, 4)


def test_encipher_kid_rsa():
    assert encipher_kid_rsa(1, (5, 2)) == 2
    assert encipher_kid_rsa(1, (8, 3)) == 3
    assert encipher_kid_rsa(1, (7, 2)) == 2


def test_decipher_kid_rsa():
    assert decipher_kid_rsa(2, (5, 3)) == 1
    assert decipher_kid_rsa(3, (8, 3)) == 1
    assert decipher_kid_rsa(2, (7, 4)) == 1


def test_encode_morse():
    assert encode_morse('ABC') == '.-|-...|-.-.'
    assert encode_morse('SMS ') == '...|--|...||'
    assert encode_morse('SMS\n') == '...|--|...||'
    assert encode_morse('') == ''
    assert encode_morse(' ') == '||'
    assert encode_morse(' ', sep='`') == '``'
    assert encode_morse(' ', sep='``') == '````'
    assert encode_morse('!@#$%^&*()_+') == '-.-.--|.--.-.|...-..-|-.--.|-.--.-|..--.-|.-.-.'


def test_decode_morse():
    assert decode_morse('-.-|.|-.--') == 'KEY'
    assert decode_morse('.-.|..-|-.||') == 'RUN'
    raises(KeyError, lambda: decode_morse('.....----'))


def test_lfsr_sequence():
    raises(TypeError, lambda: lfsr_sequence(1, [1], 1))
    raises(TypeError, lambda: lfsr_sequence([1], 1, 1))
    F = FF(2)
    assert lfsr_sequence([F(1)], [F(1)], 2) == [F(1), F(1)]
    assert lfsr_sequence([F(0)], [F(1)], 2) == [F(1), F(0)]
    F = FF(3)
    assert lfsr_sequence([F(1)], [F(1)], 2) == [F(1), F(1)]
    assert lfsr_sequence([F(0)], [F(2)], 2) == [F(2), F(0)]
    assert lfsr_sequence([F(1)], [F(2)], 2) == [F(2), F(2)]


def test_lfsr_autocorrelation():
    raises(TypeError, lambda: lfsr_autocorrelation(1, 2, 3))
    F = FF(2)
    s = lfsr_sequence([F(1), F(0)], [F(0), F(1)], 5)
    assert lfsr_autocorrelation(s, 2, 0) == 1
    assert lfsr_autocorrelation(s, 2, 1) == -1


def test_lfsr_connection_polynomial():
    F = FF(2)
    x = symbols("x")
    s = lfsr_sequence([F(1), F(0)], [F(0), F(1)], 5)
    assert lfsr_connection_polynomial(s) == x**2 + 1
    s = lfsr_sequence([F(1), F(1)], [F(0), F(1)], 5)
    assert lfsr_connection_polynomial(s) == x**2 + x + 1


def test_elgamal_private_key():
    a, b, _ = elgamal_private_key(digit=100)
    assert isprime(a)
    assert is_primitive_root(b, a)
    assert len(bin(a)) >= 102


def test_elgamal():
    dk = elgamal_private_key(5)
    ek = elgamal_public_key(dk)
    P = ek[0]
    assert P - 1 == decipher_elgamal(encipher_elgamal(P - 1, ek), dk)
    raises(ValueError, lambda: encipher_elgamal(P, dk))
    raises(ValueError, lambda: encipher_elgamal(-1, dk))


def test_dh_private_key():
    p, g, _ = dh_private_key(digit = 100)
    assert isprime(p)
    assert is_primitive_root(g, p)
    assert len(bin(p)) >= 102


def test_dh_public_key():
    p1, g1, a = dh_private_key(digit = 100)
    p2, g2, ga = dh_public_key((p1, g1, a))
    assert p1 == p2
    assert g1 == g2
    assert ga == pow(g1, a, p1)


def test_dh_shared_key():
    prk = dh_private_key(digit = 100)
    p, _, ga = dh_public_key(prk)
    b = randrange(2, p)
    sk = dh_shared_key((p, _, ga), b)
    assert sk == pow(ga, b, p)
    raises(ValueError, lambda: dh_shared_key((1031, 14, 565), 2000))


def test_padded_key():
    assert padded_key('b', 'ab') == 'ba'
    raises(ValueError, lambda: padded_key('ab', 'ace'))
    raises(ValueError, lambda: padded_key('ab', 'abba'))


def test_bifid():
    raises(ValueError, lambda: encipher_bifid('abc', 'b', 'abcde'))
    assert encipher_bifid('abc', 'b', 'abcd') == 'bdb'
    raises(ValueError, lambda: decipher_bifid('bdb', 'b', 'abcde'))
    assert encipher_bifid('bdb', 'b', 'abcd') == 'abc'
    raises(ValueError, lambda: bifid_square('abcde'))
    assert bifid5_square("B") == \
        bifid5_square('BACDEFGHIKLMNOPQRSTUVWXYZ')
    assert bifid6_square('B0') == \
        bifid6_square('B0ACDEFGHIJKLMNOPQRSTUVWXYZ123456789')


def test_playfair_matrix():
    assert playfair_matrix("PLAYFAIR") == "PLAYFIRBCDEGHKMNOQSTUVWXZ"
    assert playfair_matrix("ENCRYPT") == "ENCRYPTABDFGHIKLMOQSUVWXZ"
    assert playfair_matrix("JANGLE") == "IANGLEBCDFHKMOPQRSTUVWXYZ"


def test_encipher_playfair():
    assert encipher_playfair("ABC","ABC") == "BCHC"
    assert encipher_playfair("THEQUICKBROWNFOXJUMPEDOVERTHELAZYDOG","PLAYFAIR") == "QMHNPEKSCBQVTPSVEPEFMIVLGIQMGPFWFCVO"


def test_decipher_playfair():
    assert decipher_playfair("BCHC","ABC") == "ABCX"
    assert decipher_playfair("QMHNPEKSCBQVTPSVEPEFMIVLGIQMGPFWFCVO","PLAYFAIR") == "THEQUICKBROWNFOXIUMPEDOVERTHELAZYDOG"


def test_encipher_rijndael():
    hmsg = "6bc1bee22e409f96e93d7e117393172aae2d8a571e03ac9c9eb76fac45af8e5130c81c46a35ce411e5fbc1191a0a52eff69f2445df4f9b17ad2b417be66c3710"
    key128 = "2b7e151628aed2a6abf7158809cf4f3c"
    key192 = "8e73b0f7da0e6452c810f32b809079e562f8ead2522c6b7b"
    key256 = "603deb1015ca71be2b73aef0857d77811f352c073b6108d72d9810a30914dff4"
    iv = "000102030405060708090a0b0c0d0e0f"
    ctr = "f0f1f2f3f4f5f6f7f8f9fafbfcfdfeff"

    #ECB
    assert encipher_rijndael(hmsg,key128) == "3ad77bb40d7a3660a89ecaf32466ef97f5d3d58503b9699de785895a96fdbaaf43b1cd7f598ece23881b00e3ed0306887b0c785e27e8ad3f8223207104725dd4"
    assert encipher_rijndael(hmsg,key192) == "bd334f1d6e45f25ff712a214571fa5cc974104846d0ad3ad7734ecb3ecee4eefef7afd2270e2e60adce0ba2face6444e9a4b41ba738d6c72fb16691603c18e0e"
    assert encipher_rijndael(hmsg,key256) == "f3eed1bdb5d2a03c064b5a7e3db181f8591ccb10d410ed26dc5ba74a31362870b6ed21b99ca6f4f9f153e7b1beafed1d23304b7a39f9f3ff067d8d8f9e24ecc7"

    #CBC
    assert encipher_rijndael(hmsg,key128,mode='CBC',iv=iv) == "7649abac8119b246cee98e9b12e9197d5086cb9b507219ee95db113a917678b273bed6b8e3c1743b7116e69e222295163ff1caa1681fac09120eca307586e1a7"
    assert encipher_rijndael(hmsg,key192,mode='CBC',iv=iv) == "4f021db243bc633d7178183a9fa071e8b4d9ada9ad7dedf4e5e738763f69145a571b242012fb7ae07fa9baac3df102e008b0e27988598881d920a9e64f5615cd"
    assert encipher_rijndael(hmsg,key256,mode='CBC',iv=iv) == "f58c4c04d6e5f1ba779eabfb5f7bfbd69cfc4e967edb808d679f777bc6702c7d39f23369a9d9bacfa530e26304231461b2eb05e2c39be9fcda6c19078c6a9d1b"

    #CFB1
    assert encipher_rijndael(hmsg[0:4],key128,mode='CFB',iv=iv,s=1) == hexify("0110100010110011","binary")
    assert encipher_rijndael(hmsg[0:4],key192,mode='CFB',iv=iv,s=1) == hexify("1001001101011001","binary")
    assert encipher_rijndael(hmsg[0:4],key256,mode='CFB',iv=iv,s=1) == hexify("1001000000101001","binary")

    #CFB8
    assert encipher_rijndael(hmsg[0:36],key128,mode='CFB',iv=iv,s=8) == "3b79424c9c0dd436bace9e0ed4586a4f32b9"
    assert encipher_rijndael(hmsg[0:36],key192,mode='CFB',iv=iv,s=8) == "cda2521ef0a905ca44cd057cbf0d47a0678a"
    assert encipher_rijndael(hmsg[0:36],key256,mode='CFB',iv=iv,s=8) == "dc1f1a8520a64db55fcc8ac554844e889700"

    #CFB128
    assert encipher_rijndael(hmsg,key128,mode='CFB',iv=iv,s=128) == "3b3fd92eb72dad20333449f8e83cfb4ac8a64537a0b3a93fcde3cdad9f1ce58b26751f67a3cbb140b1808cf187a4f4dfc04b05357c5d1c0eeac4c66f9ff7f2e6"
    assert encipher_rijndael(hmsg,key192,mode='CFB',iv=iv,s=128) == "cdc80d6fddf18cab34c25909c99a417467ce7f7f81173621961a2b70171d3d7a2e1e8a1dd59b88b1c8e60fed1efac4c9c05f9f9ca9834fa042ae8fba584b09ff"
    assert encipher_rijndael(hmsg,key256,mode='CFB',iv=iv,s=128) == "dc7e84bfda79164b7ecd8486985d386039ffed143b28b1c832113c6331e5407bdf10132415e54b92a13ed0a8267ae2f975a385741ab9cef82031623d55b1e471"

    #OFB
    assert encipher_rijndael(hmsg,key128,mode='OFB',iv=iv) == "3b3fd92eb72dad20333449f8e83cfb4a7789508d16918f03f53c52dac54ed8259740051e9c5fecf64344f7a82260edcc304c6528f659c77866a510d9c1d6ae5e"
    assert encipher_rijndael(hmsg,key192,mode='OFB',iv=iv) == "cdc80d6fddf18cab34c25909c99a4174fcc28b8d4c63837c09e81700c11004018d9a9aeac0f6596f559c6d4daf59a5f26d9f200857ca6c3e9cac524bd9acc92a"
    assert encipher_rijndael(hmsg,key256,mode='OFB',iv=iv) == "dc7e84bfda79164b7ecd8486985d38604febdc6740d20b3ac88f6ad82a4fb08d71ab47a086e86eedf39d1c5bba97c4080126141d67f37be8538f5a8be740e484"

    #CTR
    assert encipher_rijndael(hmsg,key128,mode='CTR',ctr=ctr) == "874d6191b620e3261bef6864990db6ce9806f66b7970fdff8617187bb9fffdff5ae4df3edbd5d35e5b4f09020db03eab1e031dda2fbe03d1792170a0f3009cee"
    assert encipher_rijndael(hmsg,key192,mode='CTR',ctr=ctr) == "1abc932417521ca24f2b0459fe7e6e0b090339ec0aa6faefd5ccc2c6f4ce8e941e36b26bd1ebc670d1bd1d665620abf74f78a7f6d29809585a97daec58c6b050"
    assert encipher_rijndael(hmsg,key256,mode='CTR',ctr=ctr) == "601ec313775789a5b7a7f504bbf3d228f443e3ca4d62b59aca84e990cacaf5c52b0930daa23de94ce87017ba2d84988ddfc9c58db67aada613c2dd08457941a6"


def test_decipher_rijndael():
    hmsg = "6bc1bee22e409f96e93d7e117393172aae2d8a571e03ac9c9eb76fac45af8e5130c81c46a35ce411e5fbc1191a0a52eff69f2445df4f9b17ad2b417be66c3710"
    key128 = "2b7e151628aed2a6abf7158809cf4f3c"
    key192 = "8e73b0f7da0e6452c810f32b809079e562f8ead2522c6b7b"
    key256 = "603deb1015ca71be2b73aef0857d77811f352c073b6108d72d9810a30914dff4"
    iv = "000102030405060708090a0b0c0d0e0f"
    ctr = "f0f1f2f3f4f5f6f7f8f9fafbfcfdfeff"

    #ECB
    assert decipher_rijndael("3ad77bb40d7a3660a89ecaf32466ef97f5d3d58503b9699de785895a96fdbaaf43b1cd7f598ece23881b00e3ed0306887b0c785e27e8ad3f8223207104725dd4",key128) == hmsg
    assert decipher_rijndael("bd334f1d6e45f25ff712a214571fa5cc974104846d0ad3ad7734ecb3ecee4eefef7afd2270e2e60adce0ba2face6444e9a4b41ba738d6c72fb16691603c18e0e",key192) == hmsg
    assert decipher_rijndael("f3eed1bdb5d2a03c064b5a7e3db181f8591ccb10d410ed26dc5ba74a31362870b6ed21b99ca6f4f9f153e7b1beafed1d23304b7a39f9f3ff067d8d8f9e24ecc7",key256) == hmsg

    #CBC
    assert decipher_rijndael("7649abac8119b246cee98e9b12e9197d5086cb9b507219ee95db113a917678b273bed6b8e3c1743b7116e69e222295163ff1caa1681fac09120eca307586e1a7",key128,mode='CBC',iv=iv) == hmsg
    assert decipher_rijndael("4f021db243bc633d7178183a9fa071e8b4d9ada9ad7dedf4e5e738763f69145a571b242012fb7ae07fa9baac3df102e008b0e27988598881d920a9e64f5615cd",key192,mode='CBC',iv=iv) == hmsg
    assert decipher_rijndael("f58c4c04d6e5f1ba779eabfb5f7bfbd69cfc4e967edb808d679f777bc6702c7d39f23369a9d9bacfa530e26304231461b2eb05e2c39be9fcda6c19078c6a9d1b",key256,mode='CBC',iv=iv) == hmsg

    #CFB1
    assert decipher_rijndael(hexify("0110100010110011",'binary'),key128,mode='CFB',iv=iv,s=1) == hmsg[0:4]
    assert decipher_rijndael(hexify("1001001101011001","binary"),key192,mode='CFB',iv=iv,s=1) == hmsg[0:4]
    assert decipher_rijndael(hexify("1001000000101001","binary"),key256,mode='CFB',iv=iv,s=1) == hmsg[0:4]

    #CFB8
    assert decipher_rijndael("3b79424c9c0dd436bace9e0ed4586a4f32b9",key128,mode='CFB',iv=iv,s=8) == hmsg[0:36]
    assert decipher_rijndael("cda2521ef0a905ca44cd057cbf0d47a0678a",key192,mode='CFB',iv=iv,s=8) == hmsg[0:36]
    assert decipher_rijndael("dc1f1a8520a64db55fcc8ac554844e889700",key256,mode='CFB',iv=iv,s=8) == hmsg[0:36]

    #CFB128
    assert decipher_rijndael("3b3fd92eb72dad20333449f8e83cfb4ac8a64537a0b3a93fcde3cdad9f1ce58b26751f67a3cbb140b1808cf187a4f4dfc04b05357c5d1c0eeac4c66f9ff7f2e6",key128,mode='CFB',iv=iv,s=128) == hmsg
    assert decipher_rijndael("cdc80d6fddf18cab34c25909c99a417467ce7f7f81173621961a2b70171d3d7a2e1e8a1dd59b88b1c8e60fed1efac4c9c05f9f9ca9834fa042ae8fba584b09ff",key192,mode='CFB',iv=iv,s=128) == hmsg
    assert decipher_rijndael("dc7e84bfda79164b7ecd8486985d386039ffed143b28b1c832113c6331e5407bdf10132415e54b92a13ed0a8267ae2f975a385741ab9cef82031623d55b1e471",key256,mode='CFB',iv=iv,s=128) == hmsg

    #OFB
    assert decipher_rijndael("3b3fd92eb72dad20333449f8e83cfb4a7789508d16918f03f53c52dac54ed8259740051e9c5fecf64344f7a82260edcc304c6528f659c77866a510d9c1d6ae5e",key128,mode='OFB',iv=iv) == hmsg
    assert decipher_rijndael("cdc80d6fddf18cab34c25909c99a4174fcc28b8d4c63837c09e81700c11004018d9a9aeac0f6596f559c6d4daf59a5f26d9f200857ca6c3e9cac524bd9acc92a",key192,mode='OFB',iv=iv) == hmsg
    assert decipher_rijndael("dc7e84bfda79164b7ecd8486985d38604febdc6740d20b3ac88f6ad82a4fb08d71ab47a086e86eedf39d1c5bba97c4080126141d67f37be8538f5a8be740e484",key256,mode='OFB',iv=iv) == hmsg

    #CTR
    assert decipher_rijndael("874d6191b620e3261bef6864990db6ce9806f66b7970fdff8617187bb9fffdff5ae4df3edbd5d35e5b4f09020db03eab1e031dda2fbe03d1792170a0f3009cee",key128,mode='CTR',ctr=ctr) == hmsg
    assert decipher_rijndael("1abc932417521ca24f2b0459fe7e6e0b090339ec0aa6faefd5ccc2c6f4ce8e941e36b26bd1ebc670d1bd1d665620abf74f78a7f6d29809585a97daec58c6b050",key192,mode='CTR',ctr=ctr) == hmsg
    assert decipher_rijndael("601ec313775789a5b7a7f504bbf3d228f443e3ca4d62b59aca84e990cacaf5c52b0930daa23de94ce87017ba2d84988ddfc9c58db67aada613c2dd08457941a6",key256,mode='CTR',ctr=ctr) == hmsg
