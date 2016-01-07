# -*- coding: utf-8 -*-

"""
Classical ciphers and LFSRs
"""

from __future__ import print_function

from random import randrange
from string import whitespace, ascii_uppercase as uppercase, maketrans, printable

from sympy import nextprime
from sympy.core import Rational, S, Symbol
from sympy.core.numbers import igcdex
from sympy.core.compatibility import range
from sympy.matrices import Matrix
from sympy.ntheory import isprime, totient, primitive_root
from sympy.polys.domains import FF
from sympy.polys.polytools import gcd, Poly, invert
from sympy.utilities.iterables import flatten, uniq


def AZ(s=None):
    if not s:
        return uppercase
    t = type(s) is str
    if t:
        s = [s]
    rv = [check_and_join(i.upper().split(), uppercase, filter=True)
        for i in s]
    if t:
        return rv[0]
    return rv

bifid5 = AZ().replace('J', '')
bifid6 = AZ() + '0123456789'
bifid10 = printable

def padded_key(key, symbols, warn=False):
    """Return a string of the characters of ``symbols`` with those of ``key``
    appearing first, omitting characters in ``key`` that are not in ``symbols``
    unless ``warn`` is True (default=False).

    Examples
    ========

    >>> from sympy.crypto.crypto import AZ, padded_key
    >>> padded_key('rsa'.upper(), AZ())
    'RSABCDEFGHIJKLMNOPQTUVWXYZ'
    """
    if warn:
        extra = set(key) - set(symbols)
        if extra:
            raise ValueError('invalid characters: %s' % ''.join(sorted(extra)))
    key0 = ''.join([i for i in uniq(key) if i in symbols])
    return key0 + ''.join([i for i in symbols if i not in key0])


def check_and_join(phrase, symbols=None, filter=None):
    """
    Joins characters of `phrase` and if symbols are provided, raises
    an error if any character in the phrase is not in the symbols.

    Parameters
    ==========

    phrase:     string or list of strings to be returned as a string
    symbols:    iterable of characters allowed in phrase; if None,
                no checking is performed

    Examples
    ========

    >>> from sympy.crypto.crypto import check_and_join
    >>> check_and_join('a phrase')
    'a phrase'
    >>> check_and_join('a phrase'.upper().split())
    'APHRASE'
    >>> check_and_join('a phrase!'.upper().split(), 'ARE', filter=True)
    'ARAE'

    """
    rv = ''.join(''.join(phrase))
    if symbols is not None:
        symbols = check_and_join(symbols)
        missing = ''.join(list(sorted(set(rv) - set(symbols))))
        if missing:
            if not filter:
                raise ValueError('characters missing from symbols: "%s"' % missing)
            rv = rv.translate(None, missing)
    return rv


"""
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
"""

def cycle_list(k, n):
    """
    Returns the cyclic shift of the list range(n) by k.

    Examples
    ========

    >>> from sympy.crypto.crypto import cycle_list
    >>> cycle_list(3, 10)
    [3, 4, 5, 6, 7, 8, 9, 0, 1, 2]

    """
    return list(range(k, n)) + list(range(k))


def encipher_repeating(pt, key, symbols):
    """
    >>> from sympy.crypto.crypto import encipher_repeating as encode, AZ, uppercase
    >>> msg = 'attack at dawn'
    >>> key = 'lemon'
    >>> encode(*AZ((msg, key, uppercase)))
    'LXFOPVEFRNHR'
    """
    map = dict([(c, i) for i, c in enumerate(symbols)])
    key0 = [map[c] for c in key]
    N = len(map)
    k = len(key0)
    rv = []
    for i, m in enumerate(pt):
        rv.append(symbols[(map[m] + key0[i % k]) % N])
    return ''.join(rv)


def encipher_translate(pt, old, new):
    """
    >>> from sympy.crypto.crypto import encipher_translate as encode, AZ
    >>> msg = 'attack at dawn'
    >>> old = 'AT'
    >>> new = 'UB'
    >>> encode(*AZ((msg, old, new)))
    'UBBUCKUBDUWN'
    """
    return pt.translate(maketrans(old, new))


######## shift cipher examples ############


def encipher_shift(pt, key, symbols=None):
    """
    Performs shift cipher encryption on plaintext pt, and returns the ciphertext.

    Notes
    =====

    The shift cipher is also called the Caesar cipher, after
    Julius Caesar, who, according to Suetonius, used it with a
    shift of three to protect messages of military significance.
    Caesar's nephew Augustus reportedly used a similar cipher, but
    with a right shift of 1.


    ALGORITHM:

        INPUT:

            ``key``: an integer from 0 to 25 (the secret key)

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
    A = symbols or AZ()
    shift = len(A) - key % len(A)
    key = A[shift:] + A[:shift]
    return encipher_translate(pt, key, A)


######## affine cipher examples ############


def encipher_affine(pt, key, symbols=None):
    r"""
    Performs the affine cipher encryption on plaintext ``pt``, and returns the ciphertext.

    Encryption is based on the map `x \rightarrow ax+b` (mod `N`). where ``N`` is the
    number of characters in the alphabet. Decryption is based on
    the map `x \rightarrow cx+d` (mod `N`), where `c = a^{-1}` (mod `N`) and
    `d = -a^{-1}c` (mod `N`). (In particular, for the map to be invertible,
    we need `\mathrm{gcd}(a, N) = 1` and an error will be raised if this is
    not true.)

    Notes
    =====

    This is a straightforward generalization of the shift cipher with the
    added complexity of requiring 2 characters to be deciphered in order
    to recover the key.

    ALGORITHM:

        INPUT:

            ``a, b``: a pair integers, where ``gcd(a, 26) = 1`` (the secret key)

            ``m``: string of characters that appear in ``symbols``

            ``symbols``: valid characters for the plain text (default = uppercase letters)

        OUTPUT:

            ``c``: string of characters (the ciphertext message)

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
    >>> id = (1, 0)
    >>> encipher_affine(pt, id)
    'GONAVYBEATARMY'
    >>> rot1 = (1, 1)
    >>> pt = "GONAVYBEATARMY"
    >>> encipher_affine(pt, rot1)
    'HPOBWZCFBUBSNZ'
    >>> encipher_affine(pt, (3, 1))
    'TROBMVENBGBALV'
    >>> encipher_affine(_, (9, 17))
    'GONAVYBEATARMY'

    """
    A = check_and_join(AZ(symbols))
    N = len(A)
    a, b = key
    assert gcd(a, N) == 1
    B = ''.join([A[(a*i + b) % N] for i in range(N)])
    return encipher_translate(pt, A, B)


#################### substitution cipher ###########################


def encipher_substitution(pt, key, symbols=None):
    """
    Performs the substitution cipher encryption on plaintext ``pt``, and returns the ciphertext.

    Assumes the ``pt`` has only letters taken from ``symbols``.
    Assumes ``key`` is a permutation of the symbols. This function permutes the
    letters of the plaintext using the permutation given in ``key``.
    The description uses the inverse permutation.
    Note that if the permutation in key is order 2 (eg, a transposition) then
    the encryption permutation and the decryption permutation are the same.

    Notes
    =====

    This is a more general than the affine cipher in that the key can only
    be recovered by determining the mapping for each symbol. Though in
    practice, once a few symbols are recognized the mappings for other
    characters can be guessed quickly.

    Examples
    ========

    >>> from sympy.crypto.crypto import encipher_substitution, AZ
    >>> key = "BACDEFGHIJKLMNOPQRSTUVWXYZ"
    >>> pt = "go navy! beat army!"
    >>> encipher_substitution(AZ(pt), key)
    'GONBVYAEBTBRMY'
    >>> ct = 'GONBVYAEBTBRMY'
    >>> encipher_substitution(ct, key)
    'GONAVYBEATARMY'

    """
    A = check_and_join(AZ(symbols))
    key = check_and_join(key, A)
    pt = check_and_join(pt, A)
    return pt.translate(maketrans(A, key))


######################################################################
#################### Vigenère cipher examples ########################
######################################################################

def encipher_vigenere(pt, key, symbols=None):
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

          ``key``: a string of letters (the secret key)

          ``m``: string of letters (the plaintext message)

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

    >>> from sympy.crypto.crypto import encipher_vigenere, AZ
    >>> key = "encrypt"
    >>> pt = "meet me on monday"
    >>> encipher_vigenere(pt, key)
    'QRGKKTHRZQEBPR'

    """
    A = check_and_join(AZ(symbols))
    key = check_and_join(AZ(key), A)
    pt = check_and_join(AZ(pt), A)
    return encipher_repeating(pt, key, A)


def decipher_vigenere(ct, key, symbols=None):
    """
    Decode using the Vigenère cipher.

    Examples
    ========

    >>> from sympy.crypto.crypto import decipher_vigenere, AZ
    >>> key = "encrypt"
    >>> ct = "QRGK kt HRZQE BPR"
    >>> decipher_vigenere(*AZ((ct, key)))
    'MEETMEONMONDAY'
    """
    A = check_and_join(AZ(symbols))
    map = dict([(c, i) for i, c in enumerate(A)])
    N = len(A)   # normally, 26
    key0 = check_and_join(key, A)
    ct0 = check_and_join(ct, A)
    K = [map[c] for c in key0]
    n = len(K)
    C = [map[c] for c in ct0]
    P = [A[(-K[i % n] + c) % N] for i, c in enumerate(C)]
    return "".join(P)


#################### Hill cipher  ########################


def encipher_hill(pt, key, symbols=None, pad=None):
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

            ``pad``: character (default "Q") to use to make length of text be a mutiple of ``m``

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

    See Also
    ========
    decipher_hill

    """
    A = check_and_join(AZ(symbols))
    map = dict([(c, i) for i, c in enumerate(A)])
    pt0 = check_and_join(pt.upper(), A, filter=True)
    P = [map[c] for c in pt0]
    N = len(A)
    k = key.cols
    n = len(P)
    m, r = divmod(n, k)
    if r:
        P = P + [map.get(AZ(pad or 'Q'), 0)]*(k - r)
        m += 1
    C = [A[c % N] for j in range(m) for c in
        list(key*Matrix(k, 1, [P[i] for i in range(k*j, k*(j + 1))]))]
    return "".join(C)


def decipher_hill(ct, key, symbols=None):
    """
    Deciphering is the same as enciphering but using the inverse of the key matrix.

    Examples
    ========

    >>> from sympy.crypto.crypto import encipher_hill, decipher_hill, AZ
    >>> from sympy import Matrix

    >>> key = Matrix([[1, 2], [3, 5]])
    >>> encipher_hill("meet me on monday", key)
    'UEQDUEODOCTCWQ'
    >>> decipher_hill(_, key)
    'MEETMEONMONDAY'

    When the length of the plain text (stripped of invalid characters)
    is not a multiple of the key dimension extra characters will
    appear at the end of the enciphered and deciphered text. In order to
    decipher the text, those characters must be included in the text to
    be deciphered. In the following, the key has a dimension of 4 but the
    text is 2 short of being a multiple of 4 so two characters will be added.

    >>> key = Matrix([[1, 1, 1, 2], [0, 1, 1, 0], [2, 2, 3, 4], [1, 1, 0, 1]])
    >>> pt = "ST"
    >>> encipher_hill(pt, key)
    'HJEB'
    >>> decipher_hill(_, key)
    'STQQ'
    >>> encipher_hill(pt, key, pad="Z")
    'ISPK'
    >>> decipher_hill(_, key)
    'STZZ'

    If the last two characters of the cipher text were ignored in either case,
    the wrong plain text would be recovered:

    >>> decipher_hill("HD", key)
    'ORMV'
    >>> decipher_hill("IS", key)
    'UIKY'

    """
    A = check_and_join(AZ(symbols))
    map = dict([(c, i) for i, c in enumerate(A)])
    ct0 = AZ(ct)
    C = [map[c] for c in ct0]
    N = len(A)
    k = key.cols
    n = len(C)
    m, r = divmod(n, k)
    if r:
        C = C + [0]*(k - r)
        m += 1
    key_inv = key.inv_mod(N)
    P = [A[p % N] for j in range(m) for p in
        list(key_inv*Matrix(k, 1, [C[i] for i in range(k*j, k*(j + 1))]))]
    return "".join(P)


#################### Bifid cipher  ########################


def encipher_bifid(pt, key, symbols=None, filter=True):
    r"""
    Performs the Bifid cipher encryption on plaintext ``pt``, and returns the ciphertext.

    This is the version of the Bifid cipher that uses an `n \times n` Polybius square.

        INPUT:

            ``pt``: plaintext string

            ``key``: short string for key made of characters from `symbols`

            ``symbols``: `n \times n` characters defining the alphabet (default is string.printable)

            ``filter``: when True (default) characters in ``pt`` that are not in ``symbols`` are omitted

        OUTPUT:

            ciphertext (using Bifid5 cipher without spaces)

    See Also
    ========
    decipher_bifid, encipher_bifid5, encipher_bifid6

    """
    A = check_and_join(symbols) if symbols else bifid10
    n = len(A)**.5
    if n != int(n):
        raise ValueError('Length of alphabet (%s) is not a square number.' % len(A))
    N = int(n)
    pt0 = check_and_join(pt, A, filter=filter)
    long_key = check_and_join(''.join(uniq(key)), A, filter=filter) or A
    if len(long_key) < N**2:
      long_key = list(long_key) + [x for x in A if x not in long_key]
    # the fractionalization
    row_col = dict([(ch, divmod(i, N)) for i, ch in enumerate(long_key)])
    r, c = zip(*[row_col[x] for x in pt0])
    rc = r + c
    ch = dict([(i, ch) for ch, i in row_col.items()])
    ct = "".join((ch[i] for i in zip(rc[::2], rc[1::2])))
    return ct


def decipher_bifid(ct, key, symbols=None, filter=True):
    r"""
    Performs the Bifid cipher decryption on ciphertext ``ct``, and returns the plaintext.

    This is the version of the Bifid cipher that uses the `n \times n` Polybius square.

        INPUT:

            ``ct``: ciphertext string

            ``key``: short string for key made of characters from ``symbols``

            ``symbols``: `n \times n` characters defining the alphabet (default=string.printable, a `10 \times 10` matrix)

            ``filter``: when True (default) characters in ``pt`` that are not in ``symbols`` are omitted, otherwise an error is raised

        OUTPUT:

            deciphered text

    Examples
    ========

    >>> from sympy.crypto.crypto import encipher_bifid, decipher_bifid, AZ

    Do an encryption using the bifid5 alphabet:

    >>> alp = AZ().replace('J', '')
    >>> ct = AZ("meet me on monday!")
    >>> key = AZ("gold bug")
    >>> encipher_bifid(ct, key, alp)
    'IEILHHFSTSFQYE'

    When entering the text or ciphertext, spaces are ignored so it can be
    formatted as desired. Re-entering the ciphertext from the preceding,
    putting 4 characters per line and padding with an extra J, does not
    cause problems for the deciphering:

    >>> decipher_bifid('''
    ... IEILH
    ... HFSTS
    ... FQYEJ''', key, alp)
    'MEETMEONMONDAY'

    When no alphabet is given, all 100 printable characters will be used:

    >>> encipher_bifid('hello world!', '')
    'bmtwmg-bIo*w'
    >>> decipher_bifid(_, '')
    'hello world!'

    If the key is changed, a different encryption is obtained:

    >>> key = 'gold bug'
    >>> encipher_bifid('hello world!', 'gold_bug')
    'hg2sfuei7t}w'

    And if the key used to decrypt the message is not exact, the
    original text will not be perfectly obtained:

    >>> decipher_bifid(_, 'gold pug')
    'heldo~wor6d!'

    """
    A = check_and_join(symbols) if symbols else bifid10
    n = len(A)**.5
    if n != int(n):
        raise ValueError('Length of alphabet (%s) is not a square number.' % len(A))
    N = int(n)
    ct0 = check_and_join(ct, A, filter=filter)
    long_key = check_and_join(uniq(key), A, filter=filter) or A
    if len(long_key) < N**2:
        long_key = list(long_key) + [x for x in A if x not in long_key]
    # the reverse fractionalization
    row_col = dict([(ch, divmod(i, N)) for i, ch in enumerate(long_key)])
    rc = [i for c in ct0 for i in row_col[c]]
    n = len(ct0)
    rc = zip(*(rc[:n], rc[n:]))
    ch = dict([(i, ch) for ch, i in row_col.items()])
    pt = "".join((ch[i] for i in rc))
    return pt


def bifid_square(key):
    """Return characters of ``key`` arranged in a square.

    Examples
    ========

    >>> from sympy.crypto.crypto import bifid_square, AZ, padded_key, bifid5
    >>> bifid_square(AZ().replace('J', ''))
    Matrix([
    [A, B, C, D, E],
    [F, G, H, I, K],
    [L, M, N, O, P],
    [Q, R, S, T, U],
    [V, W, X, Y, Z]])

    >>> bifid_square(padded_key('GOLD BUG', bifid5))
    Matrix([
    [G, O, L, D, B],
    [U, A, C, E, F],
    [H, I, K, M, N],
    [P, Q, R, S, T],
    [V, W, X, Y, Z]])

    See Also
    ========
    padded_key
    """
    A = ''.join(uniq(key.upper()))
    n = len(A)**.5
    if n != int(n):
        raise ValueError('Length of alphabet (%s) is not a square number.' % len(A))
    n = int(n)
    f = lambda i, j: Symbol(A[n*i + j])
    M = Matrix(n, n, f)
    return M


def encipher_bifid5(pt, key):
    r"""
    Performs the Bifid cipher encryption on plaintext ``pt``, and returns the ciphertext.

    This is the version of the Bifid cipher that uses the `5 \times 5` Polybius square.
    The letter "J" is replaced with "I" unless a ``key`` of length 25 is used.

    Notes
    =====

    The Bifid cipher was invented around 1901 by Felix Delastelle.
    It is a *fractional substitution* cipher, where letters are
    replaced by pairs of symbols from a smaller alphabet. The
    cipher uses a `5 \times 5` square filled with some ordering of the alphabet,
    except that "I" is used for "J" (this is a so-called
    Polybius square; there is a `6 \times 6` analog if you add back in "j" and also
    append onto the usual 26 letter alphabet, the digits 0, 1, ..., 9).
    According to Helen Gaines' book *Cryptanalysis*, this type of cipher
    was used in the field by the German Army during World War I.

    ALGORITHM: (5x5 case)

        INPUT:

            ``pt``: plaintext string; characters not in ``key`` are ignored

            ``key``: short string for key; duplicated characters are
            ignored and if the length is less then 25 characters, it
            will be padded with other letters from the alphabet omitting
            "J". Non-alphabetic characters are ignored.

        OUTPUT:

            ciphertext (all caps, no spaces)

        STEPS:
            1. Create the `5 \times 5` Polybius square ``S`` associated to the ``k`` as
               follows:

                a) starting top left, moving left-to-right, top-to-bottom,
                   place the letters of the key into a `5 \times 5` matrix,
                b) when finished, add the letters of the alphabet
                   not in the key until the `5 \times 5` square is filled. A
                   ``key`` 25 characters will (of course) completely
                   fill the square.

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

    >>> from sympy.crypto.crypto import encipher_bifid5 as f, decipher_bifid5 as F

    "J" will be omitted unless it is replaced with "I"

    >>> key = 'a'
    >>> F(f('JOSE', key), key)
    'OSE'
    >>> F(f('JOSE'.replace("J", "I"), key), key)
    'IOSE'

    See Also
    ========
    decipher_bifid5, encipher_bifid

    """
    pt0 = AZ(pt)
    A = key0 = AZ(''.join(uniq(key)))
    if len(key0) < 25:
      A = bifid5
      key0 = list(key0) + [x for x in A if x not in key0]
    return encipher_bifid(pt0, '', key0, filter=True)


def decipher_bifid5(ct, key):
    r"""
    Performs the Bifid cipher decryption on ciphertext ``ct``, and returns the plaintext.

    This is the version of the Bifid cipher that uses the `5 \times 5` Polybius square; the
    letter "J" is ignored unless a ``key`` of length 25 is used.

    INPUT:

        ``ct``: ciphertext string

        ``key``: short string for key; duplicated characters are
        ignored and if the length is less then 25 characters, it
        will be padded with other letters from the alphabet omitting
        "J". Non-alphabetic characters are ignored.

    OUTPUT:

        plaintext from Bifid5 cipher (all caps, no spaces)

    Examples
    ========

    >>> from sympy.crypto.crypto import encipher_bifid5, decipher_bifid5
    >>> key = "gold bug"
    >>> encipher_bifid5('meet me on friday', key)
    'IEILEHFSTSFXEE'
    >>> encipher_bifid5('meet me on monday', key)
    'IEILHHFSTSFQYE'
    >>> decipher_bifid5(_, key)
    'MEETMEONMONDAY'

    """
    ct0 = AZ(ct)
    A = key0 = AZ(''.join(uniq(key)))
    if len(key0) < 25:
        A = bifid5
        key0 = list(key0) + [x for x in A if x not in key0]
    return decipher_bifid(ct0, '', key0, filter=True)


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
    return bifid_square(padded_key(key.upper(), bifid5))


def encipher_bifid6(pt, key):
    r"""
    Performs the Bifid cipher encryption on plaintext ``pt``, and returns the ciphertext.

    This is the version of the Bifid cipher that uses the `6 \times 6` Polybius square.
    Assumes alphabet of symbols is "A", ..., "Z", "0", ..., "9" but any alphabet of 36
    single letter characters can be used.

    INPUT:

        ``pt``: plaintext string (digits okay)

        ``key``: short string for key (digits okay)

    OUTPUT:

        ciphertext from Bifid cipher (all caps, no spaces)

    See Also
    ========
    decipher_bifid6, encipher_bifid

    """
    pt0 = pt.upper()
    A = key0 = ''.join(uniq(key.upper()))
    if len(A) < 36:
        A = bifid6
        key0 = check_and_join(key0, A, filter=True)
        key0 = list(key0) + [x for x in A if x not in key0]
    return encipher_bifid(pt0, '', key0, filter=True)


def decipher_bifid6(ct, key):
    r"""
    Performs the Bifid cipher decryption on ciphertext ``ct``, and returns the plaintext.

    This is the version of the Bifid cipher that uses the `6 \times 6` Polybius square.
    Assumes alphabet of symbols is "A", ..., "Z", "0", ..., "9" but any alphabet of 36
    single letter characters can be used.

    INPUT:

        ``ct``: ciphertext string (digits okay)

        ``key``: short string for key (digits okay)

    OUTPUT:

        plaintext from Bifid cipher (all caps, no spaces)

    Examples
    ========

    >>> from sympy.crypto.crypto import encipher_bifid6, decipher_bifid6
    >>> key = "gold bug"
    >>> encipher_bifid6('meet me on monday at 8am', key)
    'KFKLJJHF5MMMKTFRGPL'
    >>> decipher_bifid6(_, key)
    'MEETMEONMONDAYAT8AM'

    """
    ct0 = ct.upper()
    A = key0 = ''.join(uniq(key.upper()))
    if len(A) < 36:
        A = bifid6
        key0 = check_and_join(key0, A, filter=True)
        key0 = list(key0) + [x for x in A if x not in key0]
    return decipher_bifid(ct0, '', key0, filter=True)


def bifid6_square(key):
    r"""
    6x6 Polybius square.

    Produces the Polybius square for the `6 \times 6` Bifid cipher.
    Assumes alphabet of symbols is "A", ..., "Z", "0", ..., "9".

    Examples
    ========

    >>> from sympy.crypto.crypto import bifid6_square
    >>> key = "gold bug"
    >>> bifid6_square(key)
    Matrix([
    [G, O, L, D, B, U],
    [A, C, E, F, H, I],
    [J, K, M, N, P, Q],
    [R, S, T, V, W, X],
    [Y, Z, 0, 1, 2, 3],
    [4, 5, 6, 7, 8, 9]])
    """
    return bifid_square(padded_key(key.upper(), bifid6))

#################### RSA  #############################


def rsa_public_key(p, q, e):
    r"""
    The RSA *public key* is the pair `(n, e)`, where `n`
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

morse_char = {
    ".-": "A", "-...": "B",
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
char_morse = dict([(v, k) for k, v in morse_char.items()])


def encode_morse(pt, sep='|', mapping=None):
    """
    Encodes a plaintext into popular Morse Code with letters
    separated by `sep` and words by a double `sep`.

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

    mapping = mapping or char_morse
    assert sep not in mapping
    word_sep = 2*sep
    mapping[" "] = word_sep
    suffix = pt and pt[-1] in whitespace

    # normalize whitespace
    pt = (' ' if word_sep else '').join(pt.split())
    # omit unmapped chars
    chars = set(''.join(pt.split()))
    ok = set(mapping.keys())
    pt = pt.translate(None, ''.join(chars - ok))

    morsestring = []
    words = pt.split()
    for word in words:
        morseword = []
        for letter in word:
            morseletter = mapping[letter]
            morseword.append(morseletter)

        word = sep.join(morseword)
        morsestring.append(word)

    return word_sep.join(morsestring) + (word_sep if suffix else '')


def decode_morse(mc, sep='|', mapping=None):
    """
    Decodes a Morse Code with letters separated by `sep`
    (default is '|') and words by `word_sep` (default is '||)
    into plaintext.

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

    mapping = mapping or morse_char
    word_sep = 2*sep
    characterstring = []
    words = mc.strip(word_sep).split(word_sep)
    for word in words:
        letters = word.split(sep)
        chars = [mapping[c] for c in letters]
        word = "".join(chars)
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
