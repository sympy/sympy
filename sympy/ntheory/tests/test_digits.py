import os
import random
import re
import sys

from collections import Counter

from sympy.ntheory import (
    count_digits,
    is_palindromic_integer,
)


random.seed(os.urandom(10 ** 6))
py_version = sys.version_info


def test_count_digits():
    for i in range(10):
        base = random.choice(range(2, 101))
        m = random.choice(range(1, 101))
        digits = (
            random.choices(range(base), k=m) if py_version[0] == 3 and py_version[1] >= 6
            else [random.choice(range(base)) for i in range(m)]
        )
        while digits[0] == 0:
            random.shuffle(digits)
        N = sum(d * base ** i for i, d in enumerate(reversed(digits)))
        counter = Counter(digits)
        res_counter = count_digits(N, base=base)
        assert counter == res_counter


def test_is_palindromic_integer__base_smaller_than_2__raises_value_error():
    def base_smaller_than_2__raises_value_error():
        try:
            is_palindromic_integer(random.choice(range(10 ** 12 + 1)), random.choice(range(-(10 ** 12 + 1), 2)))
        except ValueError:
            return True
        return False
    assert base_smaller_than_2__raises_value_error()


def test_is_palindromic_integer__valid_palindrome__returns_true():
    sample = random.sample(range(10 ** 12 + 1), 10)
    for N in sample:
        st = str(N)
        pal = st + st[::-1]
        if random.choice([True, False]):
            ch = pal[int(len(pal) / 2)]
            pal = re.sub(r'{}+'.format(ch), ch, pal)
        N = int(pal)
        assert is_palindromic_integer(N) is True


def test_is_palindromic_integer__non_palindrome__returns_false():
    sample = random.sample(range(10 ** 12 + 1), 10)
    for N in sample:
        st = str(N)
        nonpal = list(st + st)
        random.shuffle(nonpal)
        while nonpal[0] == nonpal[-1]:
            random.shuffle(nonpal)
        N = int(''.join(nonpal))
        assert is_palindromic_integer(N) is False
