import random
import string
from collections import Counter

from sympy import Array
from sympy.external import import_module
from sympy.tensor.array.expressions.from_to_einsum import einsum_to_sympy_array

numpy = import_module("numpy")


def _compare_einsum_numpy(path, *a):
    sa = [Array(i) for i in a]
    assert (numpy.einsum(path, *a) == numpy.array(einsum_to_sympy_array(path, *sa).as_explicit().tolist())).all()


def _get_random_path(*a):
    dim = sum([len(i.shape) for i in a])
    letters = string.ascii_letters[:dim]
    indices1 = numpy.random.choice(list(letters), size=dim, replace=True)
    indices1_by_arg = []
    counter = 0
    for i in a:
        indices1_by_arg.append(list(indices1[counter:(counter+len(i.shape))]))
        counter += len(i.shape)
    indices_counted = Counter(indices1)
    indices_sing = [k for k, v in indices_counted.items() if v == 1]
    indices_mult = [k for k, v in indices_counted.items() if v > 1]
    min0 = 0 if indices_sing else 1
    indices2 = indices_sing + list(numpy.random.choice(indices_mult, size=numpy.random.randint(min0, len(indices_mult)), replace=False))
    random.shuffle(indices2)
    path = f"{','.join([''.join(i) for i in indices1_by_arg])}->{''.join(indices2)}"
    return path


def test_from_to_einsum_numpy():
    if numpy is None:
        return

    a = numpy.arange(3**4).reshape(3, 3, 3, 3)
    b = numpy.arange(3**4).reshape(3, 3, 3, 3)
    b *= 11
    c = numpy.arange(3**2).reshape(3, 3)

    _compare_einsum_numpy("ijkl->jkil", a)
    _compare_einsum_numpy("iijk->kj", a)
    _compare_einsum_numpy("iijk->kji", a)

    for i in range(5):
        path = _get_random_path(a, b)
        _compare_einsum_numpy(path, a, b)

    for i in range(2):
        path = _get_random_path(c, c, c)
        _compare_einsum_numpy(path, 7*c, 11*c, 13*c)
