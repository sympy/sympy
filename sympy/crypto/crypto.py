# -*- coding: utf-8 -*-

"""
Classical ciphers and LFSRs
"""

from __future__ import print_function

from random import randrange

from sympy import nextprime
from sympy.core import Rational, S, Symbol
from sympy.core.numbers import igcdex
from sympy.core.compatibility import range
from sympy.matrices import Matrix
from sympy.ntheory import isprime, totient, primitive_root
from sympy.polys.domains import FF
from sympy.polys.polytools import gcd, Poly, invert
from sympy.utilities.iterables import flatten, uniq


def alphabet_of_cipher(symbols="ABCDEFGHIJKLMNOPQRSTUVWXYZ"):
    """
    Returns the list of characters in the string input defining the alphabet.

    Notes
    =====

    First, some basic definitions.

    A *substitution cipher* is a method of encryption by which
    "units" (not necessarily characters) of plaintext are replaced with
    ciphertext according to a regular system. The "units" may be
    characters (ie, words of length `1`), words of length `2`, and so forth.

    A *transposition cipher* is a method of encryption by which
    the positions held by "units" of plaintext are replaced by a
    permutation of the plaintext. That is, the order of the units is
    changed using a bijective function on the characters' positions
    to perform the encryption.

    A *monoalphabetic cipher* uses fixed substitution over the entire
    message, whereas a *polyalphabetic cipher* uses a number of substitutions
    at different times in the message.

    Each of these ciphers require an alphabet for the messages to be
    constructed from.

    Examples
    ========

    >>> from sympy.crypto.crypto import alphabet_of_cipher
    >>> alphabet_of_cipher()
    ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
     'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
    >>> L = [str(i) for i in range(10)] + ['a', 'b', 'c']; L
    ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c']
    >>> A = "".join(L); A
    '0123456789abc'
    >>> alphabet_of_cipher(A)
    ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c']
    >>> alphabet_of_cipher()
    ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
     'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']

    """
    symbols = "".join(symbols)
    return list(symbols)


######## shift cipher examples ############


def cycle_list(k, n):
    """
    Returns the cyclic shift of the list range(n) by k.

    Examples
    ========

    >>> from sympy.crypto.crypto import cycle_list, alphabet_of_cipher
    >>> L = cycle_list(3,26); L
    [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
    17, 18, 19, 20, 21, 22, 23, 24, 25, 0, 1, 2]
    >>> A = alphabet_of_cipher(); A
    ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
     'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
    >>> [A[i] for i in L]
    ['D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
     'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'A', 'B', 'C']

    """
    L = list(range(n))
    return L[k:] + L[:k]


def encipher_shift(pt, key, symbols="ABCDEFGHIJKLMNOPQRSTUVWXYZ"):
    """
    Performs shift cipher encryption on plaintext pt, and returns the ciphertext.

    Notes
    =====

    The shift cipher is also called the Caesar cipher, after
    Julius Caesar, who, according to Suetonius, used it with a
    shift of three to protect messages of military significance.
    Caesar's nephew Augustus reportedtly used a similar cipher, but
    with a right shift of 1.


    ALGORITHM:

        INPUT:

            ``k``: an integer from 0 to 25 (the secret key)

            ``m``: string of upper-case letters (the plaintext message)

        OUTPUT:

            ``c``: string of upper-case letters (the ciphertext message)

        STEPS:
            0. Identify the alphabet A, ..., Z with the integers 0, ..., 25.
            1. Compute from the string ``m`` a list ``L1`` of corresponding
               integers.
            2. Compute from the list ``L1`` a new list ``L2``, given by
               adding ``(k mod 26)`` to each element in ``L1``.
            3. Compute from the list ``L2`` a string ``c`` of corresponding
               letters.

    Examples
    ========

    >>> from sympy.crypto.crypto import encipher_shift
    >>> pt = "GONAVYBEATARMY"
    >>> encipher_shift(pt, 1)
    'HPOBWZCFBUBSNZ'
    >>> encipher_shift(pt, 0)
    'GONAVYBEATARMY'
    >>> encipher_shift(pt, -1)
    'FNMZUXADZSZQLX'

    """
    symbols = "".join(symbols)
    A = alphabet_of_cipher(symbols)
    n = len(A)
    L = cycle_list(key, n)
    C = [A[(A.index(pt[i]) + key) % n] for i in range(len(pt))]
    return "".join(C)


######## affine cipher examples ############


def encipher_affine(pt, key, symbols="ABCDEFGHIJKLMNOPQRSTUVWXYZ"):
    r"""
    Performs the affine cipher encryption on plaintext ``pt``, and returns the ciphertext.

    Encryption is based on the map `x \rightarrow ax+b` (mod `26`). Decryption is based on
    the map `x \rightarrow cx+d` (mod `26`), where `c = a^{-1}` (mod `26`) and
    `d = -a^{-1}c` (mod `26`). (In particular, for the map to be invertible,
    we need `\mathrm{gcd}(a, 26) = 1.`)

    Notes
    =====

    This is a straightforward generalization of the shift cipher.

    ALGORITHM:

        INPUT:

            ``a, b``: a pair integers, where ``gcd(a, 26) = 1`` (the secret key)

            ``m``: string of upper-case letters (the plaintext message)

        OUTPUT:

            ``c``: string of upper-case letters (the ciphertext message)

        STEPS:
            0. Identify the alphabet "A", ..., "Z" with the integers 0, ..., 25.
            1. Compute from the string ``m`` a list ``L1`` of corresponding
               integers.
            2. Compute from the list ``L1`` a new list ``L2``, given by replacing
               ``x`` by ``a*x + b (mod 26)``, for each element ``x`` in ``L1``.
            3. Compute from the list ``L2`` a string ``c`` of corresponding
               letters.

    Examples
    ========

    >>> from sympy.crypto.crypto import encipher_affine
    >>> pt = "GONAVYBEATARMY"
    >>> encipher_affine(pt, (1, 1))
    'HPOBWZCFBUBSNZ'
    >>> encipher_affine(pt, (1, 0))
    'GONAVYBEATARMY'
    >>> pt = "GONAVYBEATARMY"
    >>> encipher_affine(pt, (3, 1))
    'TROBMVENBGBALV'
    >>> ct = "TROBMVENBGBALV"
    >>> encipher_affine(ct, (9, 17))
    'GONAVYBEATARMY'

    """
    symbols = "".join(symbols)
    A = alphabet_of_cipher(symbols)
    n = len(A)
    k1 = key[0] # multiplicative coeff "a"
    k2 = key[1] # additive coeff "b"
    L = cycle_list(k2, n)
    C = [A[(k1*A.index(pt[i]) + k2) % n] for i in range(len(pt))]
    return "".join(C)


#################### substitution cipher ###########################


def encipher_substitution(pt, key, symbols="ABCDEFGHIJKLMNOPQRSTUVWXYZ"):
    """
    Performs the substitution cipher encryption on plaintext ``pt``, and returns the ciphertext.

    Assumes the ``pt`` has only letters taken from ``symbols``.
    Assumes ``key`` is a permutation of the symbols. This function permutes the
    letters of the plaintext using the permutation given in ``key``.
    The decription uses the inverse permutation.
    Note that if the permutation in key is order 2 (eg, a transposition) then
    the encryption permutation and the decryption permutation are the same.

    Examples
    ========

    >>> from sympy.crypto.crypto import alphabet_of_cipher, encipher_substitution
    >>> symbols = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    >>> A = alphabet_of_cipher(symbols)
    >>> key = "BACDEFGHIJKLMNOPQRSTUVWXYZ"
    >>> pt = "go navy! beat army!"
    >>> encipher_substitution(pt, key)
    'GONBVYAEBTBRMY'
    >>> ct = 'GONBVYAEBTBRMY'
    >>> encipher_substitution(ct, key)
    'GONAVYBEATARMY'

    """
    symbols = "".join(symbols)
    A = alphabet_of_cipher(symbols)
    n = len(A)
    pt0 = [x.capitalize() for x in pt if x.isalnum()]
    ct = [key[A.index(x)] for x in pt0]
    return "".join(ct)


######################################################################
#################### Vigenère cipher examples ########################
######################################################################

def check_and_format(phrase, repeats=True, symbols="ABCDEFGHIJKLMNOPQRSTUVWXYZ"):
    """
    Returns a string formatted in a way as needed by Vigenère encryption.

    Parameters
    ==========

    phrase:     string to be formatted
    repeats:    boolean to allow or restrict repetition of characters
                in the formatted phrase
    symbols:    alphabet of characters allowed in phrase
                (a ValueError is returned if phrase includes
                a character that is not defined in symbols)

    Examples
    ========

    >>> from sympy.crypto.crypto import check_and_format
    >>> check_and_format('with repetitions')
    'WITHREPETITIONS'
    >>> check_and_format('without repetitions', False)
    'WITHOUREPNS'

    """
    for x in phrase:
        if x == " " or x == "\n":
            continue
        if x not in symbols.lower() + symbols.upper():
            errormsg = "Character {} of string {} is not defined in symbols {}.".format(x, phrase, symbols)
            raise ValueError(errormsg)
    rv = ''.join(''.join(phrase.split()).upper())
    if repeats:
        return rv
    return ''.join(uniq(rv))


def encipher_vigenere(pt, key, repeats = True, symbols="ABCDEFGHIJKLMNOPQRSTUVWXYZ"):
    """
    Performs the Vigenère cipher encryption on plaintext ``pt``, and returns the ciphertext.

    Notes
    =====

    The Vigenère cipher is named after Blaise de Vigenère, a sixteenth
    century diplomat and cryptographer, by a historical accident.
    Vigenère actually invented a different and more complicated cipher.
    The so-called *Vigenère cipher* was actually invented
    by Giovan Batista Belaso in 1553.

    This cipher was used in the 1800's, for example, during the American Civil War.
    The Confederacy used a brass cipher disk to implement the Vigenère cipher
    (now on display in the NSA Museum in Fort Meade) [1]_.

    The Vigenère cipher is a generalization of the shift cipher.
    Whereas the shift cipher shifts each letter by the same amount (that amount
    being the key of the shift cipher) the Vigenère cipher shifts
    a letter by an amount determined by the key (which is a word or
    phrase known only to the sender and receiver).

    For example, if the key was a single letter, such as "C", then the
    so-called Vigenere cipher is actually a shift cipher with a
    shift of `2` (since "C" is the 2nd letter of the alphabet, if
    you start counting at `0`). If the key was a word with two
    letters, such as "CA", then the so-called Vigenère cipher will
    shift letters in even positions by `2` and letters in odd positions
    are left alone (shifted by `0`, since "A" is the 0th letter, if
    you start counting at `0`).


    ALGORITHM:

        INPUT:

            ``key``: a string of upper-case letters (the secret key)

            ``m``: string of upper-case letters (the plaintext message)

            ``repeats``: True if repetition of characters in key is allowed.
                         Otherwise false.

        OUTPUT:

            ``c``: string of upper-case letters (the ciphertext message)

        STEPS:
            0. Identify the alphabet A, ..., Z with the integers 0, ..., 25.
            1. Compute from the string ``key`` a list ``L1`` of corresponding
               integers. Let ``n1 = len(L1)``.
            2. Compute from the string ``m`` a list ``L2`` of corresponding
               integers. Let ``n2 = len(L2)``.
            3. Break ``L2`` up sequencially into sublists of size ``n1``, and one sublist
               at the end of size smaller or equal to ``n1``.
            4. For each of these sublists ``L`` of ``L2``, compute a new list ``C`` given by
               ``C[i] = L[i] + L1[i] (mod 26)`` to the ``i``-th element in the sublist,
               for each ``i``.
            5. Assemble these lists ``C`` by concatenation into a new list of length ``n2``.
            6. Compute from the new list a string ``c`` of corresponding letters.

    Once it is known that the key is, say, `n` characters long, frequency analysis
    can be applied to every `n`-th letter of the ciphertext to determine the plaintext.
    This method is called *Kasiski examination* (although it was first discovered
    by Babbage).

    The cipher Vigenère actually discovered is an "auto-key" cipher
    described as follows.

    ALGORITHM:

        INPUT:

          ``key``: a string of upper-case letters (the secret key)

          ``m``: string of upper-case letters (the plaintext message)

        OUTPUT:

          ``c``: string of upper-case letters (the ciphertext message)

        STEPS:
            0. Identify the alphabet A, ..., Z with the integers 0, ..., 25.
            1. Compute from the string ``m`` a list ``L2`` of corresponding
               integers. Let ``n2 = len(L2)``.
            2. Let ``n1`` be the length of the key. Concatenate the string
               ``key`` with the first ``n2 - n1`` characters of the plaintext message.
               Compute from this string of length ``n2`` a list ``L1`` of corresponding
               integers. Note ``n2 = len(L1)``.
            3. Compute a new list ``C`` given by ``C[i] = L1[i] + L2[i] (mod 26)``,
               for each ``i``. Note ``n2 = len(C)``.
            4. Compute from the new list a string ``c`` of corresponding letters.

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Vigenere_cipher

    Examples
    ========

    Allowing character repetitions in key (default):

    >>> from sympy.crypto.crypto import encipher_vigenere
    >>> key = "encrypt"
    >>> pt = "meet me on monday"
    >>> encipher_vigenere(pt, key)
    'QRGKKTHRZQEBPR'

    Not allowing character repetitions in key:

    >>> from sympy.crypto.crypto import encipher_vigenere
    >>> key = "eeencryptt"
    >>> pt = "meet me on monday"
    >>> encipher_vigenere(pt, key, False)
    'QRGKKTHRZQEBPR'

    """
    symbols = "".join(symbols)
    A = alphabet_of_cipher(symbols)
    N = len(A)   # normally, 26
    key0 = check_and_format(key, repeats, symbols)
    K = [A.index(x) for x in key0]
    k = len(K)
    pt0 = check_and_format(pt, True, symbols)
    P = [A.index(x) for x in pt0]
    n = len(P)
    #m = n//k
    C = [(K[i % k] + P[i]) % N for i in range(n)]
    return "".join([str(A[x]) for x in C])


def decipher_vigenere(ct, key, repeats = True, symbols="ABCDEFGHIJKLMNOPQRSTUVWXYZ"):
    """
    Decode using the Vigenère cipher.

    Examples
    ========

    Allowing character repetitions in key (default):

    >>> from sympy.crypto.crypto import decipher_vigenere
    >>> key = "encrypt"
    >>> ct = "QRGK kt HRZQE BPR"
    >>> decipher_vigenere(ct, key)
    'MEETMEONMONDAY'

    Not allowing character repetitions in key:

    >>> from sympy.crypto.crypto import decipher_vigenere
    >>> key = "eeencryptt"
    >>> ct = "QRGK kt HRZQE BPR"
    >>> decipher_vigenere(ct, key, False)
    'MEETMEONMONDAY'

    """
    symbols = "".join(symbols)
    A = alphabet_of_cipher(symbols)
    N = len(A)   # normally, 26
    key0 = check_and_format(key, repeats, symbols)
    K = [A.index(x) for x in key0]
    k = len(K)
    ct0 = check_and_format(ct, True, symbols)
    C = [A.index(x) for x in ct0]
    n = len(C)
    #m = n//k
    P = [(-K[i % k] + C[i]) % N for i in range(n)]
    return "".join([str(A[x]) for x in P])


#################### Hill cipher  ########################


def encipher_hill(pt, key, symbols="ABCDEFGHIJKLMNOPQRSTUVWXYZ"):
    r"""
    Performs the Hill cipher encryption on plaintext ``pt``, and returns the ciphertext.

    Notes
    =====

    The Hill cipher [1]_, invented by Lester S. Hill in the 1920's [2]_,
    was the first polygraphic cipher in which it was practical (though barely)
    to operate on more than three symbols at once. The following discussion assumes
    an elementary knowledge of matrices.

    First, each letter is first encoded as a number. We assume here that
    "A" `\leftrightarrow` 0, "B" `\leftrightarrow` 1, ..., "Z" `\leftrightarrow` 25.
    We denote the integers `\{0, 1, ..., 25\}`
    by `Z_{26}`. Suppose your message m consists of `n` capital letters, with no spaces.
    This may be regarded an `n`-tuple M of elements of `Z_{26}`. A key in the Hill cipher
    is a `k x k` matrix `K`, all of whose entries are in `Z_{26}`, such that the matrix
    `K` is invertible (ie, that the linear transformation `K: Z_{26}^k \rightarrow Z_{26}^k`
    is one-to-one).

    ALGORITHM:

        INPUT:

            ``key``: a `k x k` invertible matrix `K`, all of whose entries are in `Z_{26}`

            ``m``: string of `n` upper-case letters (the plaintext message)
            (Note: Sage assumes that `n` is a multiple of `k`.)

        OUTPUT:

            ``c``: string of upper-case letters (the ciphertext message)

        STEPS:
            0. Identify the alphabet A, ..., Z with the integers 0, ..., 25.
            1. Compute from the string ``m`` a list ``L`` of corresponding
               integers. Let ``n = len(L)``.
            2. Break the list ``L`` up into ``t = ceiling(n/k)`` sublists
               ``L_1``, ..., ``L_t`` of size ``k`` (where the last list might be
               "padded" by 0's to ensure it is size ``k``).
            3. Compute new list ``C_1``, ..., ``C_t`` given by ``C[i] = K*L_i``
               (arithmetic is done mod 26), for each ``i``.
            4. Concatenate these into a list ``C = C_1 + ... + C_t``.
            5. Compute from ``C`` a string ``c`` of corresponding letters.
               This has length ``k*t``.

    References
    ==========

    .. [1] en.wikipedia.org/wiki/Hill_cipher
    .. [2] Lester S. Hill, Cryptography in an Algebraic Alphabet, The American
       Mathematical Monthly Vol.36, June-July 1929, pp.306-312.

    Examples
    ========

    >>> from sympy.crypto.crypto import encipher_hill
    >>> from sympy import Matrix
    >>> pt = "meet me on monday"
    >>> key = Matrix([[1, 2], [3, 5]])
    >>> encipher_hill(pt, key)
    'UEQDUEODOCTCWQ'
    >>> pt = "meet me on tuesday"
    >>> encipher_hill(pt, key)
    'UEQDUEODHBOYDJYU'
    >>> pt = "GONAVYBEATARMY"
    >>> key = Matrix([[1, 0, 1], [0, 1, 1], [2, 2, 3]])
    >>> encipher_hill(pt, key)
    'TBBYTKBEKKRLMYU'

    """
    symbols = "".join(symbols)
    A = alphabet_of_cipher(symbols)
    N = len(A)   # normally, 26
    k = key.cols
    pt0 = [x.capitalize() for x in pt if x.isalnum()]
    P = [A.index(x) for x in pt0]
    n = len(P)
    m = n//k
    if n > m*k:
        P = P + [0]*(n - m*k)
        m = m + 1
    C = [list(key*Matrix(k, 1, [P[i] for i in range(k*j, k*(j + 1))])) for j in range(m)]
    C = flatten(C)
    return "".join([A[i % N] for i in C])


def decipher_hill(ct, key, symbols="ABCDEFGHIJKLMNOPQRSTUVWXYZ"):
    """
    Deciphering is the same as enciphering but using the inverse of the key matrix.

    Examples
    ========

    >>> from sympy.crypto.crypto import decipher_hill
    >>> from sympy import Matrix
    >>> ct = "UEQDUEODOCTCWQ"
    >>> key = Matrix([[1, 2], [3, 5]])
    >>> decipher_hill(ct, key)
    'MEETMEONMONDAY'
    >>> ct = "UEQDUEODHBOYDJYU"
    >>> decipher_hill(ct, key)
    'MEETMEONTUESDAYA'

    """
    symbols = "".join(symbols)
    A = alphabet_of_cipher(symbols)
    N = len(A)   # normally, 26
    k = key.cols
    ct0 = [x.capitalize() for x in ct if x.isalnum()]
    C = [A.index(x) for x in ct0]
    n = len(C)
    m = n//k
    if n > m*k:
        C = C + [0]*(n - m*k)
        m = m + 1
    key_inv = key.inv_mod(N)
    P = [list(key_inv*Matrix(k, 1, [C[i] for i in range(k*j, k*(j + 1))])) for j in range(m)]
    P = flatten(P)
    return "".join([A[i % N] for i in P])


#################### Bifid cipher  ########################


def encipher_bifid5(pt, key):
    r"""
    Performs the Bifid cipher encryption on plaintext ``pt``, and returns the ciphertext.

    This is the version of the Bifid cipher that uses the `5 \times 5` Polybius square.

    Notes
    =====

    The Bifid cipher was invented around 1901 by Felix Delastelle.
    It is a *fractional substitution* cipher, where letters are
    replaced by pairs of symbols from a smaller alphabet. The
    cipher uses a `5 \times 5` square filled with some ordering of the alphabet,
    except that "i"s and "j"s are identified (this is a so-called
    Polybius square; there is a `6 \times 6` analog if you add back in "j" and also
    append onto the usual 26 letter alphabet, the digits 0, 1, ..., 9).
    According to Helen Gaines' book *Cryptanalysis*, this type of cipher
    was used in the field by the German Army during World War I.

    ALGORITHM: (5x5 case)

        INPUT:

            ``pt``: plaintext string (no "j"s)

            ``key``: short string for key (no repetitions, no "j"s)

        OUTPUT:

            ciphertext (using Bifid5 cipher in all caps, no spaces, no "J"s)

        STEPS:
            1. Create the `5 \times 5` Polybius square ``S`` associated to the ``k`` as
               follows:

                a) starting top left, moving left-to-right, top-to-bottom,
                   place the letters of the key into a 5x5 matrix,
                b) when finished, add the letters of the alphabet
                   not in the key until the 5x5 square is filled

            2. Create a list ``P`` of pairs of numbers which are the coordinates
               in the Polybius square of the letters in ``pt``.
            3. Let ``L1`` be the list of all first coordinates of ``P`` (length
               of ``L1 = n``), let ``L2`` be the list of all second coordinates
               of ``P`` (so the length of ``L2`` is also ``n``).
            4. Let ``L`` be the concatenation of ``L1`` and ``L2`` (length ``L = 2*n``),
               except that consecutive numbers are paired ``(L[2*i], L[2*i + 1])``.
               You can regard ``L`` as a list of pairs of length ``n``.
            5. Let ``C`` be the list of all letters which are of the form
               ``S[i, j]``, for all ``(i, j)`` in ``L``. As a string, this
               is the ciphertext ``ct``.

    Examples
    ========

    >>> from sympy.crypto.crypto import encipher_bifid5
    >>> pt = "meet me on monday"
    >>> key = "encrypt"
    >>> encipher_bifid5(pt, key)
    'LNLLQNPPNPGADK'
    >>> pt = "meet me on friday"
    >>> encipher_bifid5(pt, key)
    'LNLLFGPPNPGRSK'

    """
    A = alphabet_of_cipher()
    # first make sure the letters are capitalized
    # and text has no spaces
    key = uniq(key)
    key0 = [x.capitalize() for x in key if x.isalnum()]
    pt0 = [x.capitalize() for x in pt if x.isalnum()]
    # create long key
    long_key = key0 + [x for x in A if (not(x in key0) and x != "J")]
    n = len(pt0)
    # the fractionalization
    pairs = [[long_key.index(x)//5, long_key.index(x) % 5] for x in pt0]
    tmp_cipher = flatten([x[0] for x in pairs] + [x[1] for x in pairs])
    ct = "".join([long_key[5*tmp_cipher[2*i] + tmp_cipher[2*i + 1]] for i in range(n)])
    return ct


def decipher_bifid5(ct, key):
    r"""
    Performs the Bifid cipher decryption on ciphertext ``ct``, and returns the plaintext.

    This is the version of the Bifid cipher that uses the `5 \times 5` Polybius square.

    INPUT:

        ``ct``: ciphertext string (digits okay)

        ``key``: short string for key (no repetitions, digits okay)

    OUTPUT:

        plaintext from Bifid5 cipher (all caps, no spaces, no "J"s)

    Examples
    ========

    >>> from sympy.crypto.crypto import encipher_bifid5, decipher_bifid5
    >>> key = "encrypt"
    >>> pt = "meet me on monday"
    >>> encipher_bifid5(pt, key)
    'LNLLQNPPNPGADK'
    >>> ct = 'LNLLQNPPNPGADK'
    >>> decipher_bifid5(ct, key)
    'MEETMEONMONDAY'

    """
    A = alphabet_of_cipher()
    # first make sure the letters are capitalized
    # and text has no spaces
    key = uniq(key)
    key0 = [x.capitalize() for x in key if x.isalnum()]
    ct0 = [x.capitalize() for x in ct if x.isalnum()]
    # create long key
    long_key = key0 + [x for x in A if (not(x in key0) and x != "J")]
    n = len(ct0)
    # the fractionalization
    pairs = flatten([[long_key.index(x)//5, long_key.index(x) % 5] for x in ct0 if x != "J"])
    tmp_plain = flatten([[pairs[i], pairs[n + i]] for i in range(n)])
    pt = "".join([long_key[5*tmp_plain[2*i] + tmp_plain[2*i + 1]] for i in range(n)])
    return pt


def bifid5_square(key):
    r"""
    5x5 Polybius square.

    Produce the Polybius square for the `5 \times 5` Bifid cipher.

    Examples
    ========

    >>> from sympy.crypto.crypto import bifid5_square
    >>> bifid5_square("gold bug")
    Matrix([
    [G, O, L, D, B],
    [U, A, C, E, F],
    [H, I, K, M, N],
    [P, Q, R, S, T],
    [V, W, X, Y, Z]])

    """
    A = alphabet_of_cipher()
    # first make sure the letters are capitalized
    # and key has no spaces or duplicates
    key = uniq(key)
    key0 = [x.capitalize() for x in key if x.isalnum()]
    # create long key
    long_key = key0 + [x for x in A if (not(x in key0) and x != "J")]
    f = lambda i, j: Symbol(long_key[5*i + j])
    M = Matrix(5, 5, f)
    return M


def encipher_bifid6(pt, key):
    r"""
    Performs the Bifid cipher encryption on plaintext ``pt``, and returns the ciphertext.

    This is the version of the Bifid cipher that uses the `6 \times 6` Polybius square.
    Assumes alphabet of symbols is "A", ..., "Z", "0", ..., "9".

    INPUT:

        ``pt``: plaintext string (digits okay)

        ``key``: short string for key (no repetitions, digits okay)

    OUTPUT:

        ciphertext from Bifid cipher (all caps, no spaces)

    Examples
    ========

    >>> from sympy.crypto.crypto import encipher_bifid6
    >>> key = "encrypt"
    >>> pt = "meet me on monday at 8am"
    >>> encipher_bifid6(pt, key)
    'HNHOKNTA5MEPEGNQZYG'
    >>> encipher_bifid6(pt, key)
    'HNHOKNTA5MEPEGNQZYG'

    """
    A = alphabet_of_cipher() + [str(a) for a in range(10)]
    # first make sure the letters are capitalized
    # and text has no spaces
    key = uniq(key)
    key0 = [x.capitalize() for x in key if x.isalnum()]
    pt0 = [x.capitalize() for x in pt if x.isalnum()]
    # create long key
    long_key = key0 + [x for x in A if not(x in key0)]
    n = len(pt0)
    # the fractionalization
    pairs = [[long_key.index(x)//6, long_key.index(x) % 6] for x in pt0]
    tmp_cipher = flatten([x[0] for x in pairs] + [x[1] for x in pairs])
    ct = "".join([long_key[6*tmp_cipher[2*i] + tmp_cipher[2*i + 1]] for i in range(n)])
    return ct


def decipher_bifid6(ct, key):
    r"""
    Performs the Bifid cipher decryption on ciphertext ``ct``, and returns the plaintext.

    This is the version of the Bifid cipher that uses the `6 \times 6` Polybius square.
    Assumes alphabet of symbols is "A", ..., "Z", "0", ..., "9".

    INPUT:

        ``ct``: ciphertext string (digits okay)

        ``key``: short string for key (no repetitions, digits okay)

    OUTPUT:

        plaintext from Bifid cipher (all caps, no spaces)

    Examples
    ========

    >>> from sympy.crypto.crypto import encipher_bifid6, decipher_bifid6
    >>> key = "encrypt"
    >>> pt = "meet me on monday at 8am"
    >>> encipher_bifid6(pt, key)
    'HNHOKNTA5MEPEGNQZYG'
    >>> ct = "HNHOKNTA5MEPEGNQZYG"
    >>> decipher_bifid6(ct, key)
    'MEETMEONMONDAYAT8AM'

    """
    A = alphabet_of_cipher() + [str(a) for a in range(10)]
    # first make sure the letters are capitalized
    # and text has no spaces
    key = uniq(key)
    key0 = [x.capitalize() for x in key if x.isalnum()]
    ct0 = [x.capitalize() for x in ct if x.isalnum()]
    # create long key
    long_key = key0 + [x for x in A if not(x in key0)]
    n = len(ct0)
    # the fractionalization
    pairs = flatten([[long_key.index(x)//6, long_key.index(x) % 6] for x in ct0])
    tmp_plain = flatten([[pairs[i], pairs[n + i]] for i in range(n)])
    pt = "".join([long_key[6*tmp_plain[2*i] + tmp_plain[2*i + 1]] for i in range(n)])
    return pt


def bifid6_square(key):
    r"""
    6x6 Polybius square.

    Produces the Polybius square for the `6 \times 6` Bifid cipher.
    Assumes alphabet of symbols is "A", ..., "Z", "0", ..., "9".

    Examples
    ========

    >>> from sympy.crypto.crypto import bifid6_square
    >>> key = "encrypt"
    >>> bifid6_square(key)
    Matrix([
    [E, N, C, R, Y, P],
    [T, A, B, D, F, G],
    [H, I, J, K, L, M],
    [O, Q, S, U, V, W],
    [X, Z, 0, 1, 2, 3],
    [4, 5, 6, 7, 8, 9]])

    """
    A = alphabet_of_cipher() + [str(a) for a in range(10)]
    # first make sure the letters are capitalized
    # and text has no spaces
    key = uniq(key)
    key0 = [x.capitalize() for x in key if x.isalnum()]
    # create long key
    long_key = key0 + [x for x in A if not(x in key0)]
    f = lambda i, j: Symbol(long_key[6*i + j])
    M = Matrix(6, 6, f)
    return M


def encipher_bifid7(pt, key):
    r"""
    Performs the Bifid cipher encryption on plaintext ``pt``, and returns the ciphertext.

    This is the version of the Bifid cipher that uses the `7 \times 7` Polybius square.
    Assumes alphabet of symbols is "A", ..., "Z", "0", ..., "22".
    (Also, assumes you have some way of distinguishing "22"
    from "2", "2" juxtaposed together for deciphering...)

    INPUT:

        ``pt``: plaintext string (digits okay)

        ``key``: short string for key (no repetitions, digits okay)

    OUTPUT:

        ciphertext from Bifid7 cipher (all caps, no spaces)

    Examples
    ========

    >>> from sympy.crypto.crypto import encipher_bifid7
    >>> key = "encrypt"
    >>> pt = "meet me on monday at 8am"
    >>> encipher_bifid7(pt, key)
    'JEJJLNAA3ME19YF3J222R'

    """
    A = alphabet_of_cipher() + [str(a) for a in range(23)]
    # first make sure the letters are capitalized
    # and text has no spaces
    key = uniq(key)
    key0 = [x.capitalize() for x in key if x.isalnum()]
    pt0 = [x.capitalize() for x in pt if x.isalnum()]
    # create long key
    long_key = key0 + [x for x in A if not(x in key0)]
    n = len(pt0)
    # the fractionalization
    pairs = [[long_key.index(x)//7, long_key.index(x) % 7] for x in pt0]
    tmp_cipher = flatten([x[0] for x in pairs] + [x[1] for x in pairs])
    ct = "".join([long_key[7*tmp_cipher[2*i] + tmp_cipher[2*i + 1]] for i in range(n)])
    return ct


def bifid7_square(key):
    r"""
    7x7 Polybius square.

    Produce the Polybius square for the `7 \times 7` Bifid cipher.
    Assumes alphabet of symbols is "A", ..., "Z", "0", ..., "22".
    (Also, assumes you have some way of distinguishing "22"
    from "2", "2" juxtaposed together for deciphering...)

    Examples
    ========

    >>> from sympy.crypto.crypto import bifid7_square
    >>> bifid7_square("gold bug")
    Matrix([
    [ G,  O,  L,  D,  B,  U,  A],
    [ C,  E,  F,  H,  I,  J,  K],
    [ M,  N,  P,  Q,  R,  S,  T],
    [ V,  W,  X,  Y,  Z,  0,  1],
    [ 2,  3,  4,  5,  6,  7,  8],
    [ 9, 10, 11, 12, 13, 14, 15],
    [16, 17, 18, 19, 20, 21, 22]])

    """
    A = alphabet_of_cipher() + [str(a) for a in range(23)]
    # first make sure the letters are capitalized
    # and text has no spaces
    key = uniq(key)
    key0 = [x.capitalize() for x in key if x.isalnum()]
    # create long key
    long_key = key0 + [x for x in A if (not(x in key0))]
    f = lambda i, j: Symbol(long_key[7*i + j])
    M = Matrix(7, 7, f)
    return M


#################### RSA  #############################


def rsa_public_key(p, q, e):
    r"""
    The RSA *public key* is the pair `(n,e)`, where `n`
    is a product of two primes and `e` is relatively
    prime (coprime) to the Euler totient `\phi(n)`.

    Examples
    ========

    >>> from sympy.crypto.crypto import rsa_public_key
    >>> p, q, e = 3, 5, 7
    >>> n, e = rsa_public_key(p, q, e)
    >>> n
    15
    >>> e
    7

    """
    n = p*q
    phi = totient(n)
    if isprime(p) and isprime(q) and gcd(e, phi) == 1:
        return n, e
    return False


def rsa_private_key(p, q, e):
    r"""
    The RSA *private key* is the pair `(n,d)`, where `n`
    is a product of two primes and `d` is the inverse of
    `e` (mod `\phi(n)`).

    Examples
    ========

    >>> from sympy.crypto.crypto import rsa_private_key
    >>> p, q, e = 3, 5, 7
    >>> rsa_private_key(p, q, e)
    (15, 7)

    """
    n = p*q
    phi = totient(n)
    if isprime(p) and isprime(q) and gcd(e, phi) == 1:
        d = int(invert(e,phi))
        return n, d
    return False


def encipher_rsa(pt, puk):
    """
    In RSA, a message `m` is encrypted by computing
    `m^e` (mod `n`), where ``puk`` is the public key `(n,e)`.

    Examples
    ========

    >>> from sympy.crypto.crypto import encipher_rsa, rsa_public_key
    >>> p, q, e = 3, 5, 7
    >>> puk = rsa_public_key(p, q, e)
    >>> pt = 12
    >>> encipher_rsa(pt, puk)
    3

    """
    n, e = puk
    return pow(pt, e, n)


def decipher_rsa(ct, prk):
    """
    In RSA, a ciphertext `c` is decrypted by computing
    `c^d` (mod `n`), where ``prk`` is the private key `(n, d)`.

    Examples
    ========

    >>> from sympy.crypto.crypto import decipher_rsa, rsa_private_key
    >>> p, q, e = 3, 5, 7
    >>> prk = rsa_private_key(p, q, e)
    >>> ct = 3
    >>> decipher_rsa(ct, prk)
    12

    """
    n, d = prk
    return pow(ct, d, n)


#################### kid krypto (kid RSA) #############################


def kid_rsa_public_key(a, b, A, B):
    r"""
    Kid RSA is a version of RSA useful to teach grade school children
    since it does not involve exponentiation.

    Alice wants to talk to Bob. Bob generates keys as follows.
    Key generation:

    * Select positive integers `a, b, A, B` at random.
    * Compute `M = a b - 1`, `e = A M + a`, `d = B M + b`, `n = (e d - 1)  /M`.
    * The *public key* is `(n, e)`. Bob sends these to Alice.
    * The *private key* is `d`, which Bob keeps secret.

    Encryption: If `m` is the plaintext message then the
    ciphertext is `c = m e \pmod n`.

    Decryption: If `c` is the ciphertext message then the
    plaintext is `m = c d \pmod n`.

    Examples
    ========

    >>> from sympy.crypto.crypto import kid_rsa_public_key
    >>> a, b, A, B = 3, 4, 5, 6
    >>> kid_rsa_public_key(a, b, A, B)
    (369, 58)

    """
    M = S(a*b - 1)
    e = S(A*M + a)
    d = S(B*M + b)
    n = S((e*d - 1)//M)
    return n, e


def kid_rsa_private_key(a, b, A, B):
    """
    Compute `M = a b - 1`, `e = A M + a`, `d = B M + b`, `n = (e d - 1) / M`.
    The *private key* is `d`, which Bob keeps secret.

    Examples
    ========

    >>> from sympy.crypto.crypto import kid_rsa_private_key
    >>> a, b, A, B = 3, 4, 5, 6
    >>> kid_rsa_private_key(a, b, A, B)
    (369, 70)

    """
    M = S(a*b - 1)
    e = S(A*M + a)
    d = S(B*M + b)
    n = S((e*d - 1)//M)
    return n, d


def encipher_kid_rsa(pt, puk):
    """
    Here ``pt`` is the plaintext and ``puk`` is the public key.

    Examples
    ========

    >>> from sympy.crypto.crypto import encipher_kid_rsa, kid_rsa_public_key
    >>> pt = 200
    >>> a, b, A, B = 3, 4, 5, 6
    >>> pk = kid_rsa_public_key(a, b, A, B)
    >>> encipher_kid_rsa(pt, pk)
    161

    """
    return (pt*puk[1]) % puk[0]


def decipher_kid_rsa(ct, prk):
    """
    Here ``pt`` is the plaintext and ``prk`` is the private key.

    Examples
    ========

    >>> from sympy.crypto.crypto import kid_rsa_public_key, kid_rsa_private_key, decipher_kid_rsa, encipher_kid_rsa
    >>> a, b, A, B = 3, 4, 5, 6
    >>> d = kid_rsa_private_key(a, b, A, B)
    >>> pt = 200
    >>> pk = kid_rsa_public_key(a, b, A, B)
    >>> prk = kid_rsa_private_key(a, b, A, B)
    >>> ct = encipher_kid_rsa(pt, pk)
    >>> decipher_kid_rsa(ct, prk)
    200

    """
    n = prk[0]
    d = prk[1]
    return (ct*d) % n


#################### Morse Code ######################################


def encode_morse(pt):
    """
    Encodes a plaintext into popular Morse Code with letters separated by "|"
    and words by "||".

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Morse_code

    Examples
    ========

    >>> from sympy.crypto.crypto import encode_morse
    >>> pt = 'ATTACK THE RIGHT FLANK'
    >>> encode_morse(pt)
    '.-|-|-|.-|-.-.|-.-||-|....|.||.-.|..|--.|....|-||..-.|.-..|.-|-.|-.-'

    """

    morse_encoding_map = {"A": ".-", "B": "-...",
                     "C": "-.-.", "D": "-..",
                     "E": ".", "F": "..-.",
                     "G": "--.", "H": "....",
                     "I": "..", "J": ".---",
                     "K": "-.-", "L": ".-..",
                     "M": "--", "N": "-.",
                     "O": "---", "P": ".--.",
                     "Q": "--.-", "R": ".-.",
                     "S": "...", "T": "-",
                     "U": "..-", "V": "...-",
                     "W": ".--", "X": "-..-",
                     "Y": "-.--", "Z": "--..",
                     "0": "-----", "1": ".----",
                     "2": "..---", "3": "...--",
                     "4": "....-", "5": ".....",
                     "6": "-....", "7": "--...",
                     "8": "---..", "9": "----.",
                     ".": ".-.-.-", ",": "--..--",
                     ":": "---...", ";": "-.-.-.",
                     "?": "..--..", "-": "-...-",
                     "_": "..--.-", "(": "-.--.",
                     ")": "-.--.-", "'": ".----.",
                     "=": "-...-", "+": ".-.-.",
                     "/": "-..-.", "@": ".--.-.",
                     "$": "...-..-", "!": "-.-.--" }

    unusable_chars = "\"#%&*<>[\]^`{|}~"
    morsestring = []

    for i in unusable_chars:
        pt = pt.replace(i, "")
    pt = pt.upper()

    words = pt.split(" ")
    for word in words:
        letters = list(word)
        morseword = []
        for letter in letters:
            morseletter = morse_encoding_map[letter]
            morseword.append(morseletter)

        word = "|".join(morseword)
        morsestring.append(word)

    return "||".join(morsestring)


def decode_morse(mc):
    """
    Decodes a Morse Code with letters separated by "|"
    and words by "||" into plaintext.

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Morse_code

    Examples
    ========

    >>> from sympy.crypto.crypto import decode_morse
    >>> mc = '--|---|...-|.||.|.-|...|-'
    >>> decode_morse(mc)
    'MOVE EAST'

    """

    morse_decoding_map = {".-": "A", "-...": "B",
                          "-.-.": "C", "-..": "D",
                          ".": "E", "..-.": "F",
                          "--.": "G", "....": "H",
                          "..": "I", ".---": "J",
                          "-.-": "K", ".-..": "L",
                          "--": "M", "-.": "N",
                          "---": "O", ".--.": "P",
                          "--.-": "Q", ".-.": "R",
                          "...": "S", "-": "T",
                          "..-": "U", "...-": "V",
                          ".--": "W", "-..-": "X",
                          "-.--": "Y", "--..": "Z",
                          "-----": "0", "----": "1",
                          "..---": "2", "...--": "3",
                          "....-": "4", ".....": "5",
                          "-....": "6", "--...": "7",
                          "---..": "8", "----.": "9",
                          ".-.-.-": ".", "--..--": ",",
                          "---...": ":", "-.-.-.": ";",
                          "..--..": "?", "-...-": "-",
                          "..--.-": "_", "-.--.": "(",
                          "-.--.-": ")", ".----.": "'",
                          "-...-": "=", ".-.-.": "+",
                          "-..-.": "/", ".--.-.": "@",
                          "...-..-": "$", "-.-.--": "!"}

    characterstring = []

    if mc[-1] == "|" and mc[-2] == "|":
        mc = mc[:-2]
    words = mc.split("||")
    for word in words:
        letters = word.split("|")
        characterword = []
        for letter in letters:
            try:
                characterletter = morse_decoding_map[letter]
            except KeyError:
                return "Invalid Morse Code"
            characterword.append(characterletter)

        word = "".join(characterword)
        characterstring.append(word)
    return " ".join(characterstring)


#################### LFSRs  ##########################################


def lfsr_sequence(key, fill, n):
    r"""
    This function creates an lfsr sequence.

    INPUT:

        ``key``: a list of finite field elements,
            `[c_0, c_1, \ldots, c_k].`

        ``fill``: the list of the initial terms of the lfsr
            sequence, `[x_0, x_1, \ldots, x_k].`

        ``n``: number of terms of the sequence that the
            function returns.

    OUTPUT:

        The lfsr sequence defined by `x_{n+1} = c_k x_n + \ldots + c_0 x_{n-k}`, for
        `n \leq k`.

    Notes
    =====

    S. Golomb [G]_ gives a list of three statistical properties a
    sequence of numbers `a = \{a_n\}_{n=1}^\infty`,
    `a_n \in \{0,1\}`, should display to be considered
    "random". Define the autocorrelation of `a` to be

    .. math::

        C(k) = C(k,a) = \lim_{N\rightarrow \infty} {1\over N}\sum_{n=1}^N (-1)^{a_n + a_{n+k}}.

    In the case where `a` is periodic with period
    `P` then this reduces to

    .. math::

        C(k) = {1\over P}\sum_{n=1}^P (-1)^{a_n + a_{n+k}}.

    Assume `a` is periodic with period `P`.

    - balance:

      .. math::

        \left|\sum_{n=1}^P(-1)^{a_n}\right| \leq 1.

    - low autocorrelation:

       .. math::

         C(k) = \left\{ \begin{array}{cc} 1,& k = 0,\\ \epsilon, & k \ne 0. \end{array} \right.

      (For sequences satisfying these first two properties, it is known
      that `\epsilon = -1/P` must hold.)

    - proportional runs property: In each period, half the runs have
      length `1`, one-fourth have length `2`, etc.
      Moreover, there are as many runs of `1`'s as there are of
      `0`'s.

    References
    ==========

    .. [G] Solomon Golomb, Shift register sequences, Aegean Park Press, Laguna Hills, Ca, 1967

    Examples
    ========

    >>> from sympy.crypto.crypto import lfsr_sequence
    >>> from sympy.polys.domains import FF
    >>> F = FF(2)
    >>> fill = [F(1), F(1), F(0), F(1)]
    >>> key = [F(1), F(0), F(0), F(1)]
    >>> lfsr_sequence(key, fill, 10)
    [1 mod 2, 1 mod 2, 0 mod 2, 1 mod 2, 0 mod 2, 1 mod 2, 1 mod 2, 0 mod 2, 0 mod 2, 1 mod 2]

    """
    if not isinstance(key, list):
        raise TypeError("key must be a list")
    if not isinstance(fill, list):
        raise TypeError("fill must be a list")
    p = key[0].mod
    F = FF(p)
    s = fill
    k = len(fill)
    L = []
    for i in range(n):
        s0 = s[:]
        L.append(s[0])
        s = s[1:k]
        x = sum([int(key[i]*s0[i]) for i in range(k)])
        s.append(F(x))
    return L       # use [x.to_int() for x in L] for int version


def lfsr_autocorrelation(L, P, k):
    """
    This function computes the lsfr autocorrelation function.

    INPUT:

        ``L``: is a periodic sequence of elements of `GF(2)`.
        ``L`` must have length larger than ``P``.

        ``P``: the period of ``L``

        ``k``: an integer (`0 < k < p`)

    OUTPUT:

        the ``k``-th value of the autocorrelation of the LFSR ``L``

    Examples
    ========

    >>> from sympy.crypto.crypto import lfsr_sequence, lfsr_autocorrelation
    >>> from sympy.polys.domains import FF
    >>> F = FF(2)
    >>> fill = [F(1), F(1), F(0), F(1)]
    >>> key = [F(1), F(0), F(0), F(1)]
    >>> s = lfsr_sequence(key, fill, 20)
    >>> lfsr_autocorrelation(s, 15, 7)
    -1/15
    >>> lfsr_autocorrelation(s, 15, 0)
    1

    """
    if not isinstance(L, list):
        raise TypeError("L (=%s) must be a list" % L)
    P = int(P)
    k = int(k)
    L0 = L[:P]     # slices makes a copy
    L1 = L0 + L0[:k]
    L2 = [(-1)**(L1[i].to_int() + L1[i + k].to_int()) for i in range(P)]
    tot = sum(L2)
    return Rational(tot, P)


def lfsr_connection_polynomial(s):
    """
    This function computes the lsfr connection polynomial.

    INPUT:

        ``s``: a sequence of elements of even length, with entries in a finite field

    OUTPUT:

        ``C(x)``: the connection polynomial of a minimal LFSR yielding ``s``.

    This implements the algorithm in section 3 of J. L. Massey's article [M]_.

    References
    ==========

    .. [M] James L. Massey, "Shift-Register Synthesis and BCH Decoding."
        IEEE Trans. on Information Theory, vol. 15(1), pp. 122-127, Jan 1969.

    Examples
    ========

    >>> from sympy.crypto.crypto import lfsr_sequence, lfsr_connection_polynomial
    >>> from sympy.polys.domains import FF
    >>> F = FF(2)
    >>> fill = [F(1), F(1), F(0), F(1)]
    >>> key = [F(1), F(0), F(0), F(1)]
    >>> s = lfsr_sequence(key, fill, 20)
    >>> lfsr_connection_polynomial(s)
    x**4 + x + 1
    >>> fill = [F(1), F(0), F(0), F(1)]
    >>> key = [F(1), F(1), F(0), F(1)]
    >>> s = lfsr_sequence(key, fill, 20)
    >>> lfsr_connection_polynomial(s)
    x**3 + 1
    >>> fill = [F(1), F(0), F(1)]
    >>> key = [F(1), F(1), F(0)]
    >>> s = lfsr_sequence(key, fill, 20)
    >>> lfsr_connection_polynomial(s)
    x**3 + x**2 + 1
    >>> fill = [F(1), F(0), F(1)]
    >>> key = [F(1), F(0), F(1)]
    >>> s = lfsr_sequence(key, fill, 20)
    >>> lfsr_connection_polynomial(s)
    x**3 + x + 1

    """
    # Initialization:
    p = s[0].mod
    F = FF(p)
    x = Symbol("x")
    C = 1*x**0
    B = 1*x**0
    m = 1
    b = 1*x**0
    L = 0
    N = 0
    while N < len(s):
        if L > 0:
            dC = Poly(C).degree()
            r = min(L + 1, dC + 1)
            coeffsC = [C.subs(x, 0)] + [C.coeff(x**i) for i in range(1, dC + 1)]
            d = (s[N].to_int() + sum([coeffsC[i]*s[N - i].to_int() for i in range(1, r)])) % p
        if L == 0:
            d = s[N].to_int()*x**0
        if d == 0:
            m += 1
            N += 1
        if d > 0:
            if 2*L > N:
                C = (C - d*((b**(p - 2)) % p)*x**m*B).expand()
                m += 1
                N += 1
            else:
                T = C
                C = (C - d*((b**(p - 2)) % p)*x**m*B).expand()
                L = N + 1 - L
                m = 1
                b = d
                B = T
                N += 1
    dC = Poly(C).degree()
    coeffsC = [C.subs(x, 0)] + [C.coeff(x**i) for i in range(1, dC + 1)]
    return sum([coeffsC[i] % p*x**i for i in range(dC + 1) if coeffsC[i] is not None])


#################### ElGamal  #############################


def elgamal_private_key(digit=10):
    """
    Return three number tuple as private key.

    Elgamal encryption is based on the mathmatical problem
    called the Discrete Logarithm Problem (DLP). For example,

    `a^{b} \equiv c \pmod p`

    In general, if a and b are known, c is easily
    calculated. If b is unknown, it is hard to use
    a and c to get b.

    Parameters
    ==========

    digit : Key length in binary

    Returns
    =======

    (p, r, d) : p = prime number, r = primitive root, d = random number


    Examples
    ========

    >>> from sympy.crypto.crypto import elgamal_private_key
    >>> from sympy.ntheory import is_primitive_root, isprime
    >>> a, b, _ = elgamal_private_key()
    >>> isprime(a)
    True
    >>> is_primitive_root(b, a)
    True

    """
    p = nextprime(2**digit)
    return p, primitive_root(p), randrange(2, p)


def elgamal_public_key(prk):
    """
    Return three number tuple as public key.

    Parameters
    ==========

    prk : Tuple (p, r, e)  generated by ``elgamal_private_key``

    Returns
    =======
    (p, r, e = r**d mod p) : d is a random number in private key.

    Examples
    ========

    >>> from sympy.crypto.crypto import elgamal_public_key
    >>> elgamal_public_key((1031, 14, 636))
    (1031, 14, 212)

    """
    return prk[0], prk[1], pow(prk[1], prk[2], prk[0])


def encipher_elgamal(m, puk):
    """
    Encrypt message with public key

    m is plain text message in int. puk is
    public key (p, r, e). In order to encrypt
    a message, a random number ``a`` between ``2`` and ``p``,
    encryped message is `c_{1}` and `c_{2}`

    `c_{1} \equiv r^{a} \pmod p`

    `c_{2} \equiv m e^{a} \pmod p`

    Parameters
    ==========

    m : int of encoded message
    puk : public key

    Returns
    =======

    (c1, c2) : Encipher into two number

    Examples
    ========

    >>> from sympy.crypto.crypto import encipher_elgamal
    >>> encipher_elgamal(100, (1031, 14, 212))     # doctest: +SKIP
    (835, 271)

    """
    if m > puk[0]:
        raise ValueError('Message {} should be less than prime {}'.format(m, puk[0]))
    r = randrange(2, puk[0])
    return pow(puk[1], r, puk[0]), m * pow(puk[2], r, puk[0]) % puk[0]


def decipher_elgamal(ct, prk):
    r"""
    Decrypt message with private key

    `ct = (c_{1}, c_{2})`

    `prk = (p, r, d)`

    According to extended Eucliden theorem,
    `u c_{1}^{d} + p n = 1`

    `u \equiv 1/{{c_{1}}^d} \pmod p`

    `u c_{2} \equiv \frac{1}{c_{1}^d} c_{2} \equiv \frac{1}{r^{ad}} c_{2} \pmod p`

    `\frac{1}{r^{ad}} m e^a \equiv \frac{1}{r^{ad}} m {r^{d a}} \equiv m \pmod p`

    Examples
    ========

    >>> from sympy.crypto.crypto import decipher_elgamal
    >>> decipher_elgamal((835, 271), (1031, 14, 636))
    100

    """
    u = igcdex(ct[0] ** prk[2], prk[0])[0]
    return u * ct[1] % prk[0]


#################### Diffie-Hellman Key Exchange  #############################

def dh_private_key(digit = 10):
    """
    Return two number tuple as private key.

    Diffie-Hellman key exchange is based on the mathematical problem
    called the Discrete Logarithm Problem (see ElGamal).

    Diffie-Hellman key exchange is divided into the following steps:

    *   Alice and Bob agree on a base that consist of a prime p and a
        primitive root of p called g
    *   Alice choses a number a and Bob choses a number b where a
        and b are random numbers with 1 < a, b < p. These are their
        private keys.
    *   Alice then publicly sends Bob `g^{a} \pmod p` while Bob sends
        Alice `g^{b} \pmod p`
    *   They both raise the received value to their secretly chose number
        (a or b) and now have both as their shared key `g^{ab} \pmod p`

    Parameters
    ==========

    digit: Key length in binary

    Returns
    =======

    (p, g, a) : p = prime number, g = primitive root of p,
                a = random number in between 2 and p - 1

    Examples
    ========

    >>> from sympy.crypto.crypto import dh_private_key
    >>> from sympy.ntheory import isprime, is_primitive_root
    >>> p, g, _ = dh_private_key()
    >>> isprime(p)
    True
    >>> is_primitive_root(g, p)
    True
    >>> p, g, _ = dh_private_key(5)
    >>> isprime(p)
    True
    >>> is_primitive_root(g, p)
    True

    """
    p = nextprime(2 ** digit)
    g = primitive_root(p)
    a = randrange(2, p)
    return p, g, a

def dh_public_key(prk):
    """
    Return two number tuple as public key.

    This is the tuple that Alice sends to Bob.

    Parameters
    ==========

    prk: Tuple (p, g, a) generated by ``dh_private_key``

    Returns
    =======

    (p, g, g^a mod p) : p, g and a as in Parameters

    Examples
    ========

    >>> from sympy.crypto.crypto import dh_private_key, dh_public_key
    >>> p, g, a = dh_private_key();
    >>> _p, _g, x = dh_public_key((p, g, a))
    >>> p == _p and g == _g
    True
    >>> x == pow(g, a, p)
    True

    """
    p, g, a = prk
    return p, g, pow(g, a, p)

def dh_shared_key(puk, b):
    """
    Return int as shared key.

    This is what Bob and Alice can both calculate using the public
    keys they received from each other and their private keys.

    Parameters
    ==========

    puk: Tuple (p, g, x) generated by ``dh_public_key``
    b: Random number in the range of 2 to p - 1
       (Chosen by second key exchange member (Bob))

    Returns
    =======

    sk: int as shared key

    Examples
    ========

    >>> from sympy.crypto.crypto import dh_private_key, dh_public_key, dh_shared_key
    >>> prk = dh_private_key();
    >>> p, g, x = dh_public_key(prk);
    >>> sk = dh_shared_key((p, g, x), 1000)
    >>> sk == pow(x, 1000, p)
    True

    """
    p, _, x = puk
    if 1 >= b or b >= p:
        raise ValueError('Value of b should be greater 1 and less than prime {}'\
                .format(p))
    return pow(x, b, p)
