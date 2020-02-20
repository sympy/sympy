import os
import random
import re

random.seed(os.urandom(10 ** 6))

from itertools import chain

from sympy.ntheory.digits import (
    is_pandigital,
    is_palindromic_integer,
)


# The first 18 zerofull pandigitals - obtained from Online
# Encyclopedia of Integer Sequences (OEIS) sequence A171102
#
# https://oeis.org/A171102/list
#
sample_zerofull_pandigitals = [
    1023456789,1023456798,1023456879,1023456897,
    1023456978,1023456987,1023457689,1023457698,
    1023457869,1023457896,1023457968,1023457986,
    1023458679,1023458697,1023458769,1023458796,
    1023458967,1023458976
]
smallest_zerofull_pandigital = min(sample_zerofull_pandigitals)
sample_zeroless_pandigitals = [
    int(str(n).replace('0', ''))
    for n in sample_zerofull_pandigitals
]
smallest_zeroless_pandigital = min(sample_zeroless_pandigitals)


def test_is_pandigital__zerofull_of_insufficient_length__returns_false():
    sample_size = random.choice(range(100, 10 ** 6 + 1))
    sample = random.sample(range(smallest_zerofull_pandigital), sample_size)
    for n in sample:
        assert is_pandigital(n) is False


def test_is_pandigital__zerofull_of_insufficient_digit_frequency__returns_false():
    sample = sample_zerofull_pandigitals[::]
    ch, idx, fr = random.choice(str(1234567890)), random.choice(range(10)), random.choice(range(2, 11))
    for n in sample:
        digits = list(str(n))
        [digits.insert(idx, ch) for j in range(fr)]
        _n = int(''.join(digits))
        assert is_pandigital(_n, freq='{}+'.format(fr)) is False


def test_is_pandigital__zerofull_of_inexact_digit_frequency__returns_false():
    sample = sample_zerofull_pandigitals[::]
    ch, idx, fr = random.choice(str(1234567890)), random.choice(range(10)), random.choice(range(2, 11))
    for n in sample:
        digits = list(str(n))
        [digits.insert(idx, ch) for j in range(fr)]
        _n = int(''.join(digits))
        assert is_pandigital(_n, freq='{}+'.format(fr)) is False


def test_is_pandigital__valid_zerofull__returns_true():
    for n in sample_zerofull_pandigitals:
        assert is_pandigital(n) is True


def test_is_pandigital__valid_zerofull_of_sufficient_digit_frequency__returns_true():
    sample = sample_zerofull_pandigitals[::]
    ch, fr = random.choice(str(1234567890)), random.choice(range(2, 11))
    for n in sample:
        digits = [_d for _d in chain(*[list(str(n)) for i in range(fr)])]
        random.shuffle(digits)
        digits.insert(0, (ch if ch != '0' else '1'))
        _n = int(''.join(digits))
        assert is_pandigital(_n, freq='{}+'.format(fr)) is True


def test_is_pandigital__valid_zerofull_of_exact_digit_frequency__returns_true():
    sample = sample_zerofull_pandigitals[::]
    fr = random.choice(range(2, 11))
    for n in sample:
        digits = [_d for _d in chain(*[list(str(n)) for i in range(fr)])]
        while digits[0] == '0':
            random.shuffle(digits)
        _n = int(''.join(digits))
        assert is_pandigital(_n, freq='{}'.format(fr)) is True


def test_is_pandigital__zeroless_of_insufficient_length__returns_false():
    sample_size = random.choice(range(100, 10 ** 6 + 1))
    sample = random.sample(range(smallest_zerofull_pandigital), sample_size)
    for n in sample:
        assert is_pandigital(n) is False


def test_is_pandigital__zeroless_of_insufficient_digit_frequency__returns_false():
    sample = sample_zeroless_pandigitals[::]
    ch, idx, fr = random.choice(str(123456789)), random.choice(range(9)), random.choice(range(2, 10))
    for n in sample:
        digits = list(str(n))
        [digits.insert(idx, ch) for j in range(fr)]
        _n = int(''.join(digits))
        assert is_pandigital(_n, zeroless=True, freq='{}+'.format(fr)) is False


def test_is_pandigital__zeroless_of_inexact_digit_frequency__returns_false():
    sample = sample_zeroless_pandigitals[::]
    ch, idx, fr = random.choice(str(123456789)), random.choice(range(9)), random.choice(range(2, 10))
    for n in sample:
        digits = list(str(n))
        [digits.insert(idx, ch) for j in range(fr)]
        _n = int(''.join(digits))
        assert is_pandigital(_n, zeroless=True, freq='{}+'.format(fr)) is False


def test_is_pandigital__valid_zeroless_of_sufficient_digit_frequency__returns_true():
    sample = sample_zeroless_pandigitals[::]
    ch, fr = random.choice(str(123456789)), random.choice(range(2, 10))
    for n in sample:
        digits = [_d for _d in chain(*[list(str(n)) for i in range(fr)])]
        random.shuffle(digits)
        digits.insert(0, ch)
        _n = int(''.join(digits))
        assert is_pandigital(_n, zeroless=True, freq='{}+'.format(fr)) is True


def test_is_pandigital__valid_zeroless_of_exact_digit_frequency__returns_true():
    sample = sample_zeroless_pandigitals[::]
    fr = random.choice(range(2, 10))
    for n in sample:
        digits = [_d for _d in chain(*[list(str(n)) for i in range(fr)])]
        random.shuffle(digits)
        _n = int(''.join(digits))
        assert is_pandigital(_n, zeroless=True, freq='{}'.format(fr)) is True


def test_is_pandigital__valid_zeroless__returns_true():
    for n in sample_zeroless_pandigitals:
        assert is_pandigital(n, zeroless=True) is True


def test_is_palindromic_integer__valid_palindrome__returns_true():
    sample_size = random.choice(range(100, 10 ** 6 + 1))
    sample = random.sample(range(10 ** 12 + 1), sample_size)
    for N in sample:
        st = str(N)
        pal = st + st[::-1]
        if random.choice([True, False]):
            ch = pal[int(len(pal) / 2)]
            pal = re.sub(r'{}+'.format(ch), ch, pal)
        N = int(pal)
        assert is_palindromic_integer(N) is True


def test_is_palindromic_integer__non_palindrome__returns_false():
    sample_size = random.choice(range(100, 10 ** 6 + 1))
    sample = random.sample(range(10 ** 12 + 1), sample_size)
    for N in sample:
        st = str(N)
        nonpal = list(st + st)
        random.shuffle(nonpal)
        while nonpal[0] == nonpal[-1]:
            random.shuffle(nonpal)
        N = int(''.join(nonpal))
        assert is_palindromic_integer(N) is False
