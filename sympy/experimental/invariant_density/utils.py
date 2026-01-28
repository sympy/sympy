import random
import math
import zlib


def term_size(expr):
    """Compute syntactic term size."""
    if not expr.args:
        return 1
    return 1 + sum(term_size(a) for a in expr.args)


def kolmogorov_length(obj):
    """Approximate Kolmogorov complexity via compressed string length."""
    raw = str(obj).encode("utf-8")
    return len(zlib.compress(raw))


def random_term(variables, ops, max_size):
    """Generate a random term of bounded size."""
    if max_size <= 1 or random.random() < 0.3:
        return random.choice(variables)
    op = random.choice(ops)
    left = random.randint(1, max_size - 1)
    right = max_size - left
    return op(
        random_term(variables, ops, left),
        random_term(variables, ops, right),
    )


def wilson_interval(p_hat, n, alpha=0.05):
    """Wilson score confidence interval."""
    if n == 0:
        return (0.0, 0.0)
    z = math.sqrt(2) * math.erfcinv(alpha)
    denom = 1 + z**2 / n
    center = (p_hat + z**2 / (2 * n)) / denom
    radius = (
        z
        * math.sqrt(
            (p_hat * (1 - p_hat) / n)
            + (z**2 / (4 * n**2))
        )
        / denom
    )
    return center - radius, center + radius
