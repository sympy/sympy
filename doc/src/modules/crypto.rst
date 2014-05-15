Basic Cryptography Module
=========================

Included in this module are both block ciphers and stream ciphers.

 * Shift cipher
 * Affine cipher
 * Bifid ciphers
 * Vigenere's cipher
 * substitution ciphers
 * Hill's cipher
 * RSA
 * Kid RSA
 * linear feedback shift registers (a stream cipher)
 * ElGamal encryption

.. module:: sympy.crypto.crypto

.. autofunction:: alphabet_of_cipher

.. autofunction:: cycle_list

.. autofunction:: encipher_shift

.. autofunction:: encipher_affine

.. autofunction:: encipher_substitution

.. autofunction:: encipher_vigenere

.. autofunction:: decipher_vigenere

.. autofunction:: encipher_hill

.. autofunction:: decipher_hill

.. autofunction:: encipher_bifid5

.. autofunction:: decipher_bifid5

.. autofunction:: bifid5_square

.. autofunction:: encipher_bifid6

.. autofunction:: decipher_bifid6

.. autofunction:: bifid6_square

.. autofunction:: encipher_bifid7

.. autofunction:: bifid7_square

.. autofunction:: rsa_public_key

.. autofunction:: rsa_private_key

.. autofunction:: encipher_rsa

.. autofunction:: decipher_rsa

.. autofunction:: kid_rsa_public_key

.. autofunction:: kid_rsa_private_key

.. autofunction:: encipher_kid_rsa

.. autofunction:: decipher_kid_rsa

.. autofunction:: encode_morse

.. autofunction:: decode_morse

.. autofunction:: lfsr_sequence

.. autofunction:: lfsr_autocorrelation

.. autofunction:: lfsr_connection_polynomial

.. autofunction:: elgamal_public_key

.. autofunction:: elgamal_private_key

.. autofunction:: encipher_elgamal

.. autofunction:: decipher_elgamal
