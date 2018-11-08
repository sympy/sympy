Basic Cryptography Module
=========================

.. warning::

   This module is intended for educational purposes only. Do not use the
   functions in this module for real cryptographic applications. If you wish
   to encrypt real data, we recommend using something like the `cryptography
   <https://cryptography.io/en/latest/>`_ module.

Encryption is the process of hiding a message and a cipher is a
means of doing so. Included in this module are both block and stream
ciphers:

 * Shift cipher
 * Affine cipher
 * substitution ciphers
 * Vigenere's cipher
 * Hill's cipher
 * Bifid ciphers
 * RSA
 * Kid RSA
 * linear-feedback shift registers (for stream ciphers)
 * ElGamal encryption

In a *substitution cipher* "units" (not necessarily single characters)
of plaintext are replaced with ciphertext according to a regular system.

A *transposition cipher* is a method of encryption by which
the positions held by "units" of plaintext are replaced by a
permutation of the plaintext. That is, the order of the units is
changed using a bijective function on the position of the characters
to perform the encryption.

A *monoalphabetic cipher* uses fixed substitution over the entire
message, whereas a *polyalphabetic cipher* uses a number of
substitutions at different times in the message.

.. module:: sympy.crypto.crypto

.. autofunction:: AZ

.. autofunction:: padded_key

.. autofunction:: check_and_join

.. autofunction:: cycle_list

.. autofunction:: encipher_shift

.. autofunction:: decipher_shift

.. autofunction:: encipher_affine

.. autofunction:: decipher_affine

.. autofunction:: encipher_substitution

.. autofunction:: encipher_vigenere

.. autofunction:: decipher_vigenere

.. autofunction:: encipher_hill

.. autofunction:: decipher_hill

.. autofunction:: encipher_bifid

.. autofunction:: decipher_bifid

.. autofunction:: bifid5_square

.. autofunction:: encipher_bifid5

.. autofunction:: decipher_bifid5

.. autofunction:: bifid5_square

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

.. autofunction:: encode_morse

.. autofunction:: decode_morse

.. autofunction:: lfsr_sequence

.. autofunction:: lfsr_autocorrelation

.. autofunction:: lfsr_connection_polynomial

.. autofunction:: elgamal_public_key

.. autofunction:: elgamal_private_key

.. autofunction:: encipher_elgamal

.. autofunction:: decipher_elgamal

.. autofunction:: dh_public_key

.. autofunction:: dh_private_key

.. autofunction:: dh_shared_key

.. autofunction:: encipher_elgamal

.. autofunction:: decipher_elgamal

.. autofunction:: gm_public_key

.. autofunction:: gm_private_key

.. autofunction:: encipher_gm

.. autofunction:: decipher_gm
