from sympy.core.random import (
    random_complex_number as A,
    verify_numerically as B,
    test_derivative_numerically as C,
    _randrange as D,
    _randint as E)
from sympy.utilities.decorator import deprecated

random_complex_number = deprecated(useinstead="sympy.core.random.random_complex_number",
    deprecated_since_version="1.10", issue=22433)(A)


verify_numerically = deprecated(useinstead="sympy.core.random.verify_numerically",
    deprecated_since_version="1.10", issue=22433)(B)


test_derivative_numerically = deprecated(useinstead="sympy.core.random.test_derivative_numerically",
    deprecated_since_version="1.10", issue=22433)(C)


_randrange = deprecated(useinstead="sympy.core.random._randrange",
    deprecated_since_version="1.10", issue=22433)(D)


_randint = deprecated(useinstead="sympy.core.random._randint",
    deprecated_since_version="1.10", issue=22433)(E)
