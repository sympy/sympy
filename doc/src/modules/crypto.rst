Basic Cryptography Module
============================

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

.. module:: sympy.crypto

.. autofunction:: alphabet_of_cipher

.. autofunction:: cycle_list

.. autofunction:: encipher_shift

.. autofunction:: encipher_affine

.. autofunction:: encipher_substitution

.. autofunction:: encipher_vigenere

.. autofunction:: decipher_vigenere

.. autofunction:: matrix_inverse_mod

.. autofunction:: encipher_hill

.. autofunction:: decipher_hill

.. autofunction:: encipher_bifid5

.. autofunction:: bifid5_square

.. autofunction:: decipher_bifid5

.. autofunction:: encipher_bifid6

.. autofunction:: decipher_bifid6

.. autofunction:: bifid6_square

.. autofunction:: rsa_public_key

.. autofunction:: rsa_private_key

.. autofunction:: encipher_rsa

.. autofunction:: decipher_rsa

.. autofunction:: kid_rsa_public_key

.. autofunction:: kid_rsa_private_key

.. autofunction:: encipher_kid_rsa

.. autofunction:: decipher_kid_rsa

.. autofunction:: lfsr_sequence

.. autofunction:: lfsr_autocorrelation

.. autofunction:: lfsr_connection_polynomial
