from functools import wraps

from sympy.tensor.array.expressions import from_array_to_indexed
from sympy.utilities.exceptions import sympy_deprecation_warning


def _conv_to_from_warning():
    sympy_deprecation_warning(
        "module has been renamed by replacing 'conv_' with 'from_' in its name",
        deprecated_since_version="1.11",
        active_deprecations_target="deprecated-conv-array-expr-module-names",
    )


def _conv_to_from_decorator(f):
    @wraps(f)
    def new_f(*args, **kwargs):
        _conv_to_from_warning()
        return f(*args, **kwargs)
    return new_f


convert_array_to_indexed = _conv_to_from_decorator(from_array_to_indexed.convert_array_to_indexed)
